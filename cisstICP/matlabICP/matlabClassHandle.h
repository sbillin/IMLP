#ifndef _matlabClassHandle_h
#define _matlabClassHandle_h
#include "mex.h"
#include "stdint.h"
#include <string>
#include <cstring>
#include <typeinfo>

#include "matlabExtras.h"

#define CLASS_HANDLE_SIGNATURE 0xFF00F0A5
template<class base> class class_handle
{
public:
  class_handle(base* ptr) : ptr_m(ptr), name_m(typeid(base).name()) { signature_m = CLASS_HANDLE_SIGNATURE; }
  ~class_handle() { signature_m = 0; delete ptr_m; }
  bool isValid() { return ((signature_m == CLASS_HANDLE_SIGNATURE) && !strcmp(name_m.c_str(), typeid(base).name())); }
  base* ptr() { return ptr_m; }

private:
public:
  uint32_t signature_m;
  std::string name_m;
  base* ptr_m;
};

template<class base> inline mxArray* convertPtr2Mat(base* ptr)
{
  mexLock();
  mxArray* out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((uint64_t*)mxGetData(out)) = reinterpret_cast<uint64_t>(new class_handle<base>(ptr));
  return out;
}

template<class base> inline class_handle<base> *convertMat2HandlePtr(const mxArray *in)
{
  if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
    mexErrMsgTxt("Matlab class handle must be a real uint64 scalar.");
  class_handle<base>* ptr = reinterpret_cast<class_handle<base>*>(*((uint64_t*)mxGetData(in)));
  //std::stringstream ss;
  //ss << "ptr: " << ptr << std::endl
  //  << "signature_m: " << ptr->signature_m << std::endl
  //  << "name_m: " << ptr->name_m << std::endl
  //  << "CLASS_HANDLE_SIGNATURE: " << CLASS_HANDLE_SIGNATURE << std::endl
  //  << "signature_m == CLASS_HANDLE_SIGNATURE: " << (ptr->signature_m == CLASS_HANDLE_SIGNATURE) << std::endl
  //  << "typeid(base).name(): " << "\"" << typeid(base).name() << "\"" << std::endl
  //  << "name_m.c_str()): " << "\"" << ptr->name_m.c_str() << "\"" << std::endl
  //  << "!strcmp(name_m.c_str(), typeid(base).name()): " << !strcmp(ptr->name_m.c_str(), typeid(base).name()) << std::endl;
  //MEX_PRINT(ss.str().c_str());
  if (!ptr->isValid())
    mexErrMsgTxt("Handle not valid.");
  return ptr;
}

template<class base> inline base *convertMat2Ptr(const mxArray *in)
{
  return convertMat2HandlePtr<base>(in)->ptr();
}

template<class base> inline void destroyObject(const mxArray *in)
{
  delete convertMat2HandlePtr<base>(in);
  mexUnlock();
}

#endif // _matlabClassHandle_h

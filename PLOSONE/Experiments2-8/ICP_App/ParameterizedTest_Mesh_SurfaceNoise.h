// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// ****************************************************************************
#ifndef _ParameterizedTest_Mesh_SurfaceNoise_H
#define _ParameterizedTest_Mesh_SurfaceNoise_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include <windows.h>

#include <cisstOSAbstraction.h>

#include "cisstMesh.h"
#include "cisstCovTree_Mesh.h"
#include "cisstAlgorithmICP_StdICP_Mesh.h"
#include "cisstAlgorithmICP_RobustICP_Mesh.h"
#include "cisstAlgorithmICP_IMLP_Mesh.h"
#include "cisstAlgorithmICP_IMLP_MahalDist_Mesh.h"
#include "cisstAlgorithmICP_IMLP_ClosestPoint_Mesh.h"

#include "ParameterizedTest_PointCloud_SurfaceNoise.h"
#include "ParameterizedTest.h"
#include "utility.h"

// declerations: test routines
void Run_ParameterizedTest_Mesh_SurfaceNoise(TestParameters params);
void Run_ParameterizedTest_Mesh_SurfaceNoise_Outliers(TestParameters params);

#endif
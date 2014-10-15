//
//   Copyright 2014 DreamWorks Animation LLC.
//
//   Licensed under the Apache License, Version 2.0 (the "Apache License")
//   with the following modification; you may not use this file except in
//   compliance with the Apache License and the following modification to it:
//   Section 6. Trademarks. is deleted and replaced with:
//
//   6. Trademarks. This License does not grant permission to use the trade
//      names, trademarks, service marks, or product names of the Licensor
//      and its affiliates, except as required to comply with Section 4(c) of
//      the License and to reproduce the content of the NOTICE file.
//
//   You may obtain a copy of the Apache License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the Apache License with the above modification is
//   distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//   KIND, either express or implied. See the Apache License for the specific
//   language governing permissions and limitations under the Apache License.
//
#include "../sdc/type.h"

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Sdc {

//
//  Specializations for TypeTraits<TYPE_LOOP>:
//
template <>
Split
TypeTraits<TYPE_LOOP>::TopologicalSplitType() {
    return SPLIT_TO_TRIS;
}

template <>
int
TypeTraits<TYPE_LOOP>::LocalNeighborhoodSize() {
    return 1;
}

template <>
int
TypeTraits<TYPE_LOOP>::RegularVertexValence() {
    return 6;
}

template <>
int
TypeTraits<TYPE_LOOP>::RegularFaceValence() {
    return 3;
}

template <>
char const*
TypeTraits<TYPE_LOOP>::Label() {
    //  Might need to declare static here to keep all compilers happy...
    return "loop";
}

} // end namespace sdc

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv

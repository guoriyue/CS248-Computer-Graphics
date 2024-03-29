//
// /Users/guomingfei/Desktop/Cardinal3D-master-lab3/src/rays/ispc_bvh.h
// (Header automatically generated by the ispc compiler.)
// DO NOT EDIT THIS FILE.
//

#pragma once
#include <stdint.h>



#ifdef __cplusplus
namespace ispc { /* namespace */
#endif // __cplusplus

#ifndef __ISPC_ALIGN__
#if defined(__clang__) || !defined(_MSC_VER)
// Clang, GCC, ICC
#define __ISPC_ALIGN__(s) __attribute__((aligned(s)))
#define __ISPC_ALIGNED_STRUCT__(s) struct __ISPC_ALIGN__(s)
#else
// Visual Studio
#define __ISPC_ALIGN__(s) __declspec(align(s))
#define __ISPC_ALIGNED_STRUCT__(s) __ISPC_ALIGN__(s) struct
#endif
#endif

#ifndef __ISPC_STRUCT_Vec3__
#define __ISPC_STRUCT_Vec3__
struct Vec3 {
    float x;
    float y;
    float z;
};
#endif

#ifndef __ISPC_STRUCT_Ray__
#define __ISPC_STRUCT_Ray__
struct Ray {
    struct Vec3 point;
    struct Vec3 dir;
};
#endif

#ifndef __ISPC_STRUCT_BBox__
#define __ISPC_STRUCT_BBox__
struct BBox {
    struct Vec3 min;
    struct Vec3 max;
};
#endif

#ifndef __ISPC_STRUCT_Vec2__
#define __ISPC_STRUCT_Vec2__
struct Vec2 {
    float x;
    float y;
};
#endif


///////////////////////////////////////////////////////////////////////////
// Functions exported from ispc code
///////////////////////////////////////////////////////////////////////////
#if defined(__cplusplus) && (! defined(__ISPC_NO_EXTERN_C) || !__ISPC_NO_EXTERN_C )
extern "C" {
#endif // __cplusplus
    extern void bbox_hit(const struct Ray &ray, struct BBox * ispc_bboxs, struct Vec2 * ispc_times, bool * ispc_hits);
#if defined(__cplusplus) && (! defined(__ISPC_NO_EXTERN_C) || !__ISPC_NO_EXTERN_C )
} /* end extern C */
#endif // __cplusplus


#ifdef __cplusplus
} /* namespace */
#endif // __cplusplus

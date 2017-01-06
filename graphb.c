//#include "stdafx.h"
#ifndef __GRAPHBLAS_H__
#define __GRAPHBLAS_H__
#include <stdint.h> 

typedef int GrB_HANDLE;
typedef GrB_HANDLE GrB_Vector;
typedef GrB_HANDLE GrB_Matrix;
typedef uint32_t GrB_index;
typedef int32_t GrB_info;
typedef int32_t GrB_type;
typedef int32_t GrB_field;
typedef GrB_HANDLE GrB_Monoid;
typedef GrB_HANDLE GrB_Semiring;
typedef GrB_HANDLE GrB_Descriptor;
typedef GrB_HANDLE GrB_BinaryOp;
typedef GrB_HANDLE GrB_UnaryOp;
typedef GrB_HANDLE GrB_Variable;
#define GrB_SUCCESS 1
#define GrB_INT32 2
#define GrB_BOOL 3
#define GrB_NULL 0
#define GrB_MASK 0
#define GrB_SCMP 7
#define GrB_ALL ((GrB_index *)1)
#define GrB_NOOBJECT 8
#define GrB_PANIC 9
#define GrB_OUTOFMEM 10
#define GrB_NOMATRIX 6
#define GrB_INDEX_OUTOFBOUNDS 4
#define GrB_NOVECTOR 5
#define GrB_NOMONOID 11
#define GrB_MISMATCHEDTYPE 12
#define GrB_MISMATCHEDSIZES 13
#define GrB_ROWCOLUMN 14
#define GrB_NOUNARYOP 15
#define GrB_UNARYOPFAIL 16
#define GrB_LOR -1
#define GrB_LAND -2
#define GrB_PLUS -3
GrB_info GrB_printf(GrB_HANDLE h);
GrB_info GrB_free(GrB_HANDLE object);

GrB_info GrB_Matrix_new(GrB_Matrix *A, GrB_type d, GrB_index m, GrB_index n);
GrB_info GrB_Matrix_nrows(GrB_index *n, const GrB_Matrix A);
GrB_info GrB_Matrix_ncols(GrB_index *n, const GrB_Matrix A);
GrB_info GrB_Matrix_assign(GrB_Matrix u, int v);
GrB_info GrB_Matrix_assignij(GrB_Matrix u, int v, GrB_index i, GrB_index j);
GrB_info GrB_Matrix_nnz(GrB_index *nnz, GrB_Vector u);

GrB_info GrB_Vector_new(GrB_Vector *u, GrB_type d, GrB_index n);
GrB_info GrB_Vector_assign(GrB_Vector u, int v);
GrB_info GrB_Vector_assigni(GrB_Vector u, int v, GrB_index i);
GrB_info GrB_Vector_size(GrB_index *size, GrB_Vector u);
GrB_info GrB_Vector_nnz(GrB_index *nnz, GrB_Vector u);
GrB_info GrB_Vector_apply(GrB_Vector u, const GrB_Vector mask, const GrB_BinaryOp accum,
                          const GrB_UnaryOp op, const GrB_Vector v, 
                          const GrB_Descriptor desc, void *context);
GrB_info GrB_Vector_reduce(GrB_Variable u, const GrB_BinaryOp accum, const GrB_Monoid op,
                           const GrB_Vector A, const GrB_Descriptor desc);
GrB_info GrB_Vector_assignop(GrB_Vector u, const GrB_Vector mask, const GrB_BinaryOp accum, int v,
                             const GrB_index *indices, const GrB_index n, const GrB_Descriptor desc);
GrB_info GrB_vxm(GrB_Vector u, const GrB_Vector mask, const GrB_BinaryOp accum, const GrB_Semiring op,
                 GrB_Vector v, const GrB_Matrix A, const GrB_Descriptor desc);

GrB_info GrB_Monoid_new(GrB_Monoid *monoid, GrB_type d1, GrB_BinaryOp binary_op, int v);
GrB_info GrB_Semiring_new(GrB_Semiring  *semiring, GrB_Monoid add_op, GrB_BinaryOp mul_op);
GrB_info GrB_Descriptor_new(GrB_Descriptor *desc);
GrB_info GrB_Descriptor_set(GrB_Descriptor desc, GrB_field field, int v);
GrB_info GrB_Variable_new(GrB_Variable *c, GrB_type d, int v);

GrB_info GrB_UnaryOp_new(GrB_UnaryOp *unary_op, GrB_type d1, GrB_type d2, GrB_info (*unary_func)(GrB_Variable, const GrB_Variable, void*));
GrB_info GrB_BinaryOp_new(GrB_BinaryOp *binary_op, GrB_type d1, GrB_type d2, GrB_type d3, 
  GrB_info(*binary_func)(GrB_Variable, const GrB_Variable, const GrB_Variable, void*));

#define GrB_bool(c) (GrB_handles[c]->p.variable.v)
#define GrB_set_int32(x, c) {GrB_handles[x]->p.variable.v = c;}
#endif //GRAPHBLAS_H

//Prototype implementation of enough of GraphBLAS C99 standard to implement the first two examples
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
//#include "GraphBLAS.h"
#ifndef strcpy_s
#define strcpy_s(d, n, s) strcpy(d, s)
#endif

#define GrB_max_objects 50
#define GrB_obj_matrix 1
#define GrB_obj_vector 2
#define GrB_obj_monoid 3
#define GrB_obj_semiring 4
#define GrB_obj_descriptor 5
#define GrB_obj_variable 6
#define GrB_obj_binaryop 7
#define GrB_obj_unaryop 8

typedef struct {
  int type;
  union {
    struct {
      GrB_type type;
      GrB_index m, n;
      void *a;
    } matrix;
    struct {
      GrB_type type;
      GrB_index n;
      void *a;
    } vector;
    struct {
      GrB_type type;
      GrB_BinaryOp binary_op;
      int v;
    } monoid;
    struct {
      GrB_Monoid monoid;
      GrB_BinaryOp binary_op;
    } semiring;
    struct {
      int vals[3];
    } descriptor;
    struct {
      GrB_type type;
      int v;
    } variable;
    struct {
      GrB_type type;
      GrB_Variable result;
      GrB_Variable arg1;
      GrB_info (*f)(GrB_Variable, const GrB_Variable, void*);
    } unaryop;
    struct {
      GrB_type type;
      GrB_Variable result;
      GrB_Variable arg1;
      GrB_Variable arg2;
      GrB_info (*f)(GrB_Variable, const GrB_Variable, const GrB_Variable, void*);
    } binaryop;
  } p;
} GrB_Object;

GrB_Object *GrB_handles[GrB_max_objects];

GrB_HANDLE GrB_Handle_new() {
  int i;
  for (i = 1; i<GrB_max_objects; i++)
  if (GrB_handles[i] == NULL) {
    GrB_handles[i] = (GrB_Object *)malloc(sizeof(GrB_Object));
    return i;
  }
  return -1;
}

GrB_info GrB_free(GrB_HANDLE object) {
  GrB_Object *p;
  if (object<0 || object >= GrB_max_objects) return GrB_NOOBJECT;
  p = GrB_handles[object];
  switch (p->type) {
  case GrB_obj_variable:
  case GrB_obj_descriptor:
  case GrB_obj_monoid:
  case GrB_obj_semiring:
    break;
  case GrB_obj_matrix:
    if (p->p.matrix.a != NULL) free(p->p.matrix.a);
    break;
  case GrB_obj_vector:
    if (p->p.vector.a != NULL) free(p->p.vector.a);
    break;
  case GrB_obj_unaryop:
    GrB_free(p->p.unaryop.result);
    GrB_free(p->p.unaryop.arg1);
    break;
  case GrB_obj_binaryop:
    GrB_free(p->p.binaryop.result);
    GrB_free(p->p.binaryop.arg1);
    GrB_free(p->p.binaryop.arg2);
    break;
  default: return GrB_PANIC;
  }
  GrB_handles[object] = NULL;
  return GrB_SUCCESS;
}

int GrB_type_isvalid(GrB_type d) {
  return d == GrB_BOOL || d == GrB_INT32;
}

int GrB_type_size(GrB_type d) {
  switch (d) {
  case GrB_BOOL: return sizeof(char);
  case GrB_INT32: return sizeof(int);
  }
  return -1;
}

void GrB_type_tostring(char *s, GrB_type t) {
  switch (t) {
  case GrB_BOOL: strcpy_s(s, 10, "bool"); break;
  case GrB_INT32: strcpy_s(s, 10, "int32"); break;
  }
}

void GrB_binop_tostring(char *s, GrB_BinaryOp o) {
  if (o<0 && o>=10)
  switch (o) {
  case GrB_LOR: strcpy_s(s, 10, "lor"); break;
  case GrB_LAND: strcpy_s(s, 10, "land"); break;
  }
}

GrB_info GrB_Matrix_new(GrB_Matrix *A, GrB_type d, GrB_index m, GrB_index n) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (A == NULL || !GrB_type_isvalid(d)) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_matrix;
  p->p.matrix.m = m;
  p->p.matrix.n = n;
  p->p.matrix.type = d;
  p->p.matrix.a = malloc(GrB_type_size(d)*m*n);
  if (p->p.matrix.a == NULL) {
    GrB_free(h);
    return GrB_OUTOFMEM;
  }
  *A = h;
  return GrB_SUCCESS;
}

GrB_info GrB_Matrix_nrows(GrB_index *n, const GrB_Matrix A) {
  GrB_Object *p;
  if (n == NULL || A<0 || A >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[A];
  if (p->type != GrB_obj_matrix) return GrB_NOMATRIX;
  *n = p->p.matrix.m;
  return GrB_SUCCESS;
}
GrB_info GrB_Matrix_ncols(GrB_index *n, const GrB_Matrix A) {
  GrB_Object *p;
  if (n == NULL || A<0 || A >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[A];
  if (p->type != GrB_obj_matrix) return GrB_NOMATRIX;
  *n = p->p.matrix.n;
  return GrB_SUCCESS;
}

GrB_info GrB_printf(GrB_HANDLE h) {
  int i, j, m, n;
  signed char *pb;
  int *pi;
  GrB_Object *p;
  if (h<0 || h >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[h];
  switch (p->type) {
  case GrB_obj_matrix:
    if (p->p.matrix.a == NULL) return GrB_NOMATRIX;
    m = p->p.matrix.m;
    n = p->p.matrix.n;
    pb = (signed char *)p->p.matrix.a;
    for (i = 0; i<m; i++) {
      for (j = 0; j<n; j++)
        printf("%4d ", pb[i*n + j]);
      printf("\n");
    }
    break;
  case GrB_obj_vector:
    if (p->p.vector.a == NULL) return GrB_NOVECTOR;
    n = p->p.vector.n;
    pb = (signed char *)p->p.vector.a;
    pi = (int *)p->p.vector.a;
    for (j = 0; j < n; j++) {
      switch (p->p.vector.type) {
      case GrB_BOOL: printf("%4d ", pb[j]); break;
      case GrB_INT32: printf("%6d ", pi[j]); break;
      }
    }
    printf("\n");
    break;
  case GrB_obj_monoid:
    {
    char t[10], b[10];
    GrB_type_tostring(&t[0], p->p.monoid.type);
    GrB_binop_tostring(&b[0], p->p.monoid.binary_op);
    printf("monoid type=%s binop=%s id=%d\n", t, b, p->p.monoid.v);
    }
    break;
  case GrB_obj_semiring:
    {
    printf("semiring ");
    GrB_printf(p->p.semiring.monoid);
    GrB_printf(p->p.semiring.binary_op);
    }
    break;
  case GrB_obj_variable:
    {
    printf("%d \n", p->p.variable.v);
    }
    break;
  default: return GrB_PANIC;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Matrix_nnz(GrB_index *nnz, GrB_Matrix u) {
  int i, j, m, n;
  signed char *pb;
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  if (p->type != GrB_obj_matrix) return GrB_NOMATRIX;
  if (p->p.matrix.a == NULL) return GrB_NOMATRIX;
  m = p->p.matrix.m;
  n = p->p.matrix.n;
  pb = (signed char *)p->p.matrix.a;
  *nnz = 0;
  for (i = 0; i<m; i++) {
    for (j = 0; j<n; j++) if (pb[i*n + j]) *nnz += 1;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Matrix_assign(GrB_Matrix u, int v) {
  int i, j, m, n;
  signed char *pb;
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  switch (p->type) {
  case GrB_obj_matrix:
    if (p->p.matrix.a == NULL) return GrB_NOMATRIX;
    m = p->p.matrix.m;
    n = p->p.matrix.n;
    pb = (signed char *)p->p.matrix.a;
    for (i = 0; i<m; i++) {
      for (j = 0; j<n; j++) pb[i*n + j] = v;
    }
    break;
  default: return GrB_PANIC;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Matrix_assignij(GrB_Matrix u, int v, GrB_index i, GrB_index j) {
  GrB_index m, n;
  signed char *pb;
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  switch (p->type) {
  case GrB_obj_matrix:
    if (p->p.matrix.a == NULL) return GrB_NOMATRIX;
    m = p->p.matrix.m;
    n = p->p.matrix.n;
    pb = (signed char *)p->p.matrix.a;
    if (i >= m || j >= n) return GrB_INDEX_OUTOFBOUNDS;
    pb[i*n + j] = v;
    break;
  default: return GrB_PANIC;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Vector_new(GrB_Vector *u, GrB_type d, GrB_index n) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (u == NULL || !GrB_type_isvalid(d)) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_vector;
  p->p.vector.n = n;
  p->p.vector.type = d;
  p->p.vector.a = malloc(GrB_type_size(d)*n);
  if (p->p.vector.a == NULL) {
    GrB_free(h);
    return GrB_OUTOFMEM;
  }
  *u = h;
  return GrB_SUCCESS;
}

GrB_info GrB_Vector_assign(GrB_Vector u, int v) {
  int j, n;
  signed char *pb;
  int *pi;
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  switch (p->type) {
  case GrB_obj_vector:
    if (p->p.vector.a == NULL) return GrB_NOVECTOR;
    n = p->p.vector.n;
    pb = (signed char *)p->p.vector.a;
    pi = (int *)p->p.vector.a;
    for (j = 0; j<n; j++) {
      switch (p->p.vector.type) {
      case GrB_BOOL: pb[j] = v; break;
      case GrB_INT32: pi[j] = v; break;
      }
    }
    break;
  default: return GrB_PANIC;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Vector_size(GrB_index *size, GrB_Vector u) {
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  if (p->type != GrB_obj_vector) return GrB_NOVECTOR;
  return p->p.vector.n;
}

GrB_info GrB_Vector_nnz(GrB_index *nnz, GrB_Vector u) {
  int j;
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  if (p->type != GrB_obj_vector) return GrB_NOVECTOR;
  int n = p->p.vector.n;
  signed char *pb = (signed char *)p->p.vector.a;
  int *pi = (int *)p->p.vector.a;
  *nnz = 0;
  for (j = 0; j<n; j++) {
    switch (p->p.vector.type) {
    case GrB_BOOL: if (pb[j]) *nnz+=1; break;
    case GrB_INT32: if (pi[j]) *nnz+=1; break;
    }
  }
return GrB_SUCCESS;
}

GrB_info GrB_Vector_assigni(GrB_Vector u, int v, GrB_index i) {
  GrB_index n;
  signed char *pb;
  GrB_Object *p;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  switch (p->type) {
  case GrB_obj_vector:
    if (p->p.vector.a == NULL) return GrB_NOVECTOR;
    n = p->p.vector.n;
    pb = (signed char *)p->p.vector.a;
    if (i >= n) return GrB_INDEX_OUTOFBOUNDS;
    pb[i] = v;
    break;
  default: return GrB_PANIC;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Monoid_new(GrB_Monoid *monoid, GrB_type d, GrB_BinaryOp binary_op, int v) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (monoid == NULL || !GrB_type_isvalid(d)) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_monoid;
  p->p.monoid.type = d;
  p->p.monoid.binary_op = binary_op;
  p->p.monoid.v = v;
  *monoid = h;
  return GrB_SUCCESS;
}

GrB_info GrB_Semiring_new(GrB_Semiring  *semiring, GrB_Monoid add_op, GrB_BinaryOp mul_op) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (semiring == NULL) return GrB_PANIC;
  if (add_op<0 || add_op >= GrB_max_objects) return GrB_NOMONOID;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_semiring;
  p->p.semiring.monoid = add_op;
  p->p.semiring.binary_op = mul_op;
  *semiring = h;
  return GrB_SUCCESS;
}

GrB_info GrB_Descriptor_new(GrB_Descriptor *desc) {
  int i;
  GrB_HANDLE h;
  GrB_Object *p;
  if (desc == NULL) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_descriptor;
  for (i = 0; i<sizeof(p->p.descriptor.vals) / sizeof(int); i++) p->p.descriptor.vals[i]=0;
  *desc = h;
  return GrB_SUCCESS;
}

GrB_info GrB_Descriptor_set(GrB_Descriptor desc, GrB_field field, int v) {
  GrB_Object *p;
  if (desc<0 || desc >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[desc];
  p->p.descriptor.vals[field] = v;
  return GrB_SUCCESS;
}

GrB_info GrB_Vector_assignop(GrB_Vector u, const GrB_Vector mask, const GrB_BinaryOp accum, int v,
  const GrB_index *indices, const GrB_index n, const GrB_Descriptor desc) {
  unsigned int i;
  GrB_Object *p;
  GrB_Object *pm;
  signed char *pb;
  if (accum != GrB_NULL || desc != GrB_NULL || u == GrB_NULL || indices == NULL) return GrB_PANIC;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[u];
  if (p->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (mask<0 || mask >= GrB_max_objects) return GrB_PANIC;
  pm = GrB_handles[mask];
  if (pm->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (pm->p.vector.n != p->p.vector.n) return GrB_MISMATCHEDSIZES;
  if (p->p.vector.type != GrB_INT32) return GrB_MISMATCHEDTYPE;
  if (pm->p.vector.type != GrB_BOOL) return GrB_MISMATCHEDTYPE;
  int *pi = (int*)p->p.vector.a;
  pb = (signed char *)pm->p.vector.a;
  for (i = 0; i < p->p.vector.n; i++) {
    if (pb[i]) pi[i] = v;
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Vector_reduce(GrB_Variable c, const GrB_BinaryOp accum, const GrB_Monoid op,
  const GrB_Vector v, const GrB_Descriptor desc) {
  unsigned int i;
  GrB_Object *p;
  GrB_Object *pc;
  signed char *pb;
  if (accum != GrB_NULL || desc != GrB_NULL) return GrB_PANIC;
  if (v<0 || v >= GrB_max_objects) return GrB_PANIC;
  p = GrB_handles[v];
  if (p->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (c<0 || c >= GrB_max_objects) return GrB_PANIC;
  pc = GrB_handles[c];
  if (pc->type != GrB_obj_variable) return GrB_NOVECTOR;
  if (p->p.vector.type != pc->p.variable.type) return GrB_MISMATCHEDTYPE;
  pb = (signed char *)p->p.vector.a;
  for (i = 0; i < p->p.vector.n; i++) {
    if (*pb++) {
      pc->p.variable.v = 1;
      return GrB_SUCCESS;
    }
  } //for i
  pc->p.variable.v = 0;
  return GrB_SUCCESS;
}

GrB_info GrB_vxm(GrB_Vector u, const GrB_Vector mask, const GrB_BinaryOp accum, const GrB_Semiring op,
  GrB_Vector v, const GrB_Matrix A, const GrB_Descriptor desc) {
  unsigned int i, j, m, n;
  GrB_Object *pu;
  GrB_Object *pv;
  GrB_Object *pm;
  GrB_Object *pa;
  signed char *pU, *pV, *pA;
  int *pi;
  if (accum != GrB_NULL) return GrB_PANIC;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  pu = GrB_handles[u];
  if (pu->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (v<0 || v >= GrB_max_objects) return GrB_PANIC;
  pv = GrB_handles[v];
  if (pv->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (mask<0 || mask >= GrB_max_objects) return GrB_PANIC;
  pm = GrB_handles[mask];
  if (pm->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (A<0 || A >= GrB_max_objects) return GrB_PANIC;
  pa = GrB_handles[A];
  if (pa->type != GrB_obj_matrix) return GrB_NOMATRIX;
  if (pv->p.vector.n != pa->p.matrix.m) return GrB_ROWCOLUMN;
  if (pu->p.vector.n != pa->p.matrix.n) return GrB_MISMATCHEDSIZES;
  if (pu->p.vector.n != pm->p.vector.n) return GrB_MISMATCHEDSIZES;
  if (pu->p.vector.type != pv->p.vector.type) return GrB_MISMATCHEDTYPE;
  if (pu->p.vector.type != pa->p.vector.type) return GrB_MISMATCHEDTYPE;
  pU = (signed char *)pu->p.vector.a;
  pV = (signed char *)pv->p.vector.a;
  pi = (int *)pm->p.vector.a;
  pA = (signed char *)pa->p.matrix.a;
  m = pa->p.matrix.m;
  n = pa->p.matrix.n;
  for (j = 0; j < m; j++) {
    if (pV[j]) {
      pU[j] = 0;
      for (i = 0; i < n; i++) {
        if (pA[j*n + i] && pi[i] == 0) pU[i] = 1;
      }
    }
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Vector_apply(GrB_Vector u, const GrB_Vector mask, const GrB_BinaryOp accum,
  const GrB_UnaryOp op, const GrB_Vector v, const GrB_Descriptor desc, void *context) {
  GrB_Object *pu, *pv, *po, *pArg, *pResult;
  unsigned int i;
  signed char *pb;
  int *pi;
  if (mask != GrB_NULL || desc != GrB_NULL || accum!=GrB_PLUS) return GrB_PANIC;
  if (u<0 || u >= GrB_max_objects) return GrB_PANIC;
  pu = GrB_handles[u];
  if (pu->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (v<0 || v >= GrB_max_objects) return GrB_PANIC;
  pv = GrB_handles[v];
  if (pv->type != GrB_obj_vector) return GrB_NOVECTOR;
  if (op<0 || op >= GrB_max_objects) return GrB_PANIC;
  po = GrB_handles[op];
  if (po->type != GrB_obj_unaryop) return GrB_NOUNARYOP;
  if (pu->p.vector.type != po->p.unaryop.type) GrB_MISMATCHEDTYPE;
  if (pu->p.vector.n != pv->p.vector.n) return GrB_MISMATCHEDSIZES;
  pb = (signed char *)pv->p.vector.a;
  pi = (int *)pu->p.vector.a;
  pArg = GrB_handles[po->p.unaryop.arg1];
  pResult = GrB_handles[po->p.unaryop.result];
  for (i = 0; i < pu->p.vector.n; i++) {
    GrB_info info;
    pArg->p.variable.v = pb[i];
    info = po->p.unaryop.f(po->p.unaryop.result, po->p.unaryop.arg1, context);
    if (info != GrB_SUCCESS) return GrB_UNARYOPFAIL;
    switch (accum) {
    case GrB_PLUS : pi[i] += pResult->p.variable.v; break;
    }
  }
  return GrB_SUCCESS;
}

GrB_info GrB_Variable_new(GrB_Variable *c, GrB_type d, int v) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (c == NULL) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_variable;
  p->p.variable.type = d;
  p->p.variable.v = v;
  *c = h;
  return GrB_SUCCESS;
}

GrB_info GrB_UnaryOp_new(GrB_UnaryOp *unary_op, GrB_type d1, GrB_type d2, GrB_info (*unary_func)(GrB_Variable, const GrB_Variable, void*)) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (unary_op == NULL) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_unaryop;
  p->p.unaryop.type = d1;
  GrB_Variable_new(&p->p.unaryop.result, d1, 0);
  GrB_Variable_new(&p->p.unaryop.arg1, d2, 0);
  p->p.unaryop.f = unary_func;
  *unary_op = h;
  return GrB_SUCCESS;
}

GrB_info GrB_BinaryOp_new(GrB_BinaryOp *binary_op, GrB_type d1, GrB_type d2, GrB_type d3,
  GrB_info (*binary_func)(GrB_Variable, const GrB_Variable, const GrB_Variable, void*)) {
  GrB_HANDLE h;
  GrB_Object *p;
  if (binary_op == NULL) return GrB_PANIC;
  h = GrB_Handle_new();
  if (h < 0) return GrB_PANIC;
  p = GrB_handles[h];
  p->type = GrB_obj_binaryop;
  p->p.binaryop.type = d1;
  GrB_Variable_new(&p->p.binaryop.result, d1, 0);
  GrB_Variable_new(&p->p.binaryop.arg1, d2, 0);
  GrB_Variable_new(&p->p.binaryop.arg2, d3, 0);
  p->p.binaryop.f = binary_func;
  *binary_op = h;
  return GrB_SUCCESS;
}

/*
* Given a Boolean n x n adjacency matrix A and a source vertex s,
*  performs a BFS traversal of the graph and sets v[i] to the level
*  in which vertex i is visited (v[s] == 1).
* If i is not reachable from s, then v[i] = 0.
* (Vector v is overwritten.)
*/
//#define APPLY
#define REDUCE
GrB_info returnLevel(GrB_Variable result, const GrB_Variable arg1, void *context) {
  GrB_set_int32(result, GrB_bool(arg1) ? (int)context : 0);
  return GrB_SUCCESS;
}
GrB_info BFS(GrB_Vector *v, GrB_Matrix A, GrB_index s)  {
  GrB_index n;                       // n = # of rows of A
  GrB_Vector q;                      // vertices visited in each level
  GrB_Monoid Lor;                   // Logical or monoid
  GrB_Semiring Boolean;             // Boolean semiring <bool , bool , bool ,|| ,&& , false >
  GrB_Descriptor desc;              // Descriptor for vxm
  int32_t d = 1;                    // d = level in BFS traversal
  GrB_Variable succ;                // succ == true when some successor found
  GrB_Matrix_ncols(&n, A);
  GrB_Vector_new(v, GrB_INT32, n);     // Vector<int32 t> v(n)
  GrB_Vector_assign(*v, 0);                   // v = 0
  GrB_Vector_new(&q, GrB_BOOL, n);      // Vector<bool> q(n)
  GrB_Vector_assign(q, false);
  GrB_Vector_assigni(q, true, s);             // q[s] = true , false everywhere else
  GrB_Monoid_new(&Lor, GrB_BOOL, GrB_LOR, false);
  GrB_Semiring_new(&Boolean, Lor, GrB_LAND);
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_set(desc, GrB_MASK, GrB_SCMP);  // invert the mask

  // BFS traversal and label the vertices.
#ifdef APPLY
  {
  GrB_UnaryOp applyLevel;
  GrB_index nnz;
  int level = 0;
  GrB_UnaryOp_new(&applyLevel, GrB_INT32, GrB_BOOL, returnLevel);
  do {
      ++level; // next level(start with 1)
      GrB_Vector_apply(*v, GrB_NULL, GrB_PLUS, applyLevel, q, GrB_NULL, (void *)level); // v[q] = level
      GrB_printf(*v);
      GrB_vxm(q, *v, GrB_NULL, Boolean, q, A, desc);  // q[!v] = q ||.&& A ; finds all the
                                                      // successors from current q
      GrB_printf(q); 
      GrB_Vector_nnz(&nnz, q);
  } while (nnz); // if there is no scuccessor in q, we are done.
  GrB_free(applyLevel);
  }
#endif
#ifdef REDUCE
  GrB_Variable_new(&succ, GrB_BOOL, false);
  do {
    GrB_Vector_assignop(*v, q, GrB_NULL, d, GrB_ALL, n, GrB_NULL);  // v[q] = d
    GrB_printf(*v);
    GrB_vxm(q, *v, GrB_NULL, Boolean, q, A, desc);  // q[!v] = q ||.&& A ; finds all the
                                                    // successors from current q
    GrB_printf(q);
    GrB_Vector_reduce(succ, GrB_NULL, Lor, q, GrB_NULL); // succ = ||(q)
    GrB_printf(succ);
    d++;                                     // next level
  } while (GrB_bool(succ));                  // if there is no successor in q, we are done.
  GrB_free(succ);                            // constant no longer needed
#endif
  GrB_free(q);                               // q vector no longer needed
  GrB_free(Lor);                             // Logical or monoid no longer needed
  GrB_free(Boolean);                         // Boolean semiring no longer needed
  GrB_free(desc);                           // descriptor no longer needed
  return GrB_SUCCESS;
}

int main(void) {
  GrB_Matrix m;
  GrB_Vector v;
  GrB_Matrix_new(&m, GrB_BOOL, 14, 14);
  GrB_Matrix_assign(m, false);
  GrB_Matrix_assignij(m, true, 1, 0); //directed graph
  GrB_Matrix_assignij(m, true, 2, 0);
  GrB_Matrix_assignij(m, true, 3, 2);
  GrB_Matrix_assignij(m, true, 4, 1);
  GrB_Matrix_assignij(m, true, 5, 2);
  GrB_Matrix_assignij(m, true, 5, 4);
  GrB_Matrix_assignij(m, true, 6, 3);
  GrB_Matrix_assignij(m, true, 7, 5);
  GrB_Matrix_assignij(m, true, 9, 6);
  GrB_Matrix_assignij(m, true, 9, 8);
  GrB_Matrix_assignij(m, true, 10, 7);
  GrB_Matrix_assignij(m, true, 10, 8);
  GrB_Matrix_assignij(m, true, 11, 9);
  GrB_Matrix_assignij(m, true, 12, 10);
  GrB_Matrix_assignij(m, true, 13, 11);
  GrB_Matrix_assignij(m, true, 13, 12);
  GrB_printf(m);
  BFS(&v, m, 13);
  GrB_printf(v);
  GrB_free(m);
  GrB_free(v);
  return 0;
}

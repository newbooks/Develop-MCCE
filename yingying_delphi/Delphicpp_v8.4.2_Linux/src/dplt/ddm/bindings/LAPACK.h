#ifndef DDM__BINDINGS__LAPACK_H__
#define DDM__BINDINGS__LAPACK_H__

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus
  void dgemm_ ( void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *, void *, void *);
  void sgemm_ ( void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *, void *, void *);
  void zgemm_ ( void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *, void *, void *);
  void cgemm_ ( void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *, void *, void *);

  void dgetrf_ (void *, void *, void *, void *, void *, void *);
  void sgetrf_ (void *, void *, void *, void *, void *, void *);
  void zgetrf_ (void *, void *, void *, void *, void *, void *);
  void cgetrf_ (void *, void *, void *, void *, void *, void *);

  void dgetrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *);
  void sgetrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *);
  void zgetrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *);
  void cgetrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *);

  void dgesv_ (void *, void *, void *, void *, void *, void *, void *,
               void *);
  void sgesv_ (void *, void *, void *, void *, void *, void *, void *,
               void *);
  void zgesv_ (void *, void *, void *, void *, void *, void *, void *,
               void *);
  void cgesv_ (void *, void *, void *, void *, void *, void *, void *,
               void *);

  // Banded matrix routines
  void dgbtrf_ (void *, void *, void *, void *, void *, void *, void *,
                void *);
  void sgbtrf_ (void *, void *, void *, void *, void *, void *, void *,
                void *);
  void zgbtrf_ (void *, void *, void *, void *, void *, void *, void *,
                void *);
  void cgbtrf_ (void *, void *, void *, void *, void *, void *, void *,
                void *);

  void dgbtrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *);
  void sgbtrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *);
  void zgbtrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *);
  void cgbtrs_ (void *, void *, void *, void *, void *, void *, void *,
                void *, void *, void *, void *);

  void dgbsv_ (void *, void *, void *, void *, void *, void *, void *,
               void *, void *, void *);
  void sgbsv_ (void *, void *, void *, void *, void *, void *, void *,
               void *, void *, void *);
  void zgbsv_ (void *, void *, void *, void *, void *, void *, void *,
               void *, void *, void *);
  void cgbsv_ (void *, void *, void *, void *, void *, void *, void *,
               void *, void *, void *);
#ifdef __cplusplus
}
#endif // __cplusplus

#endif  // DDM__BINDINGS__LAPACK_H__

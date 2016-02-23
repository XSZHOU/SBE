#ifndef scattering_H
#define scattering_H

extern "C"{
/* 2d cubature called by all processes */
int IntegrandHF2d_ndiag(unsigned ndim, const double *x,  void *fdata, unsigned fdim, double *fval);
int IntegrandHF2d_diag (unsigned ndim, const double *x,  void *fdata, unsigned fdim, double *fval);
int IntegrandCC2d_diag (unsigned ndim, const double *x,  void *fdata, unsigned fdim, double *fval);
/* 1d cquad called by IntegrandCC_diag */
double Integrand1d_cquad1( const double x , void *userdata ); 
double Integrand1d_cquad2( const double x , void *userdata );
double Integrand1d_cquad3( const double x , void *userdata );
double Integrand1d_cquad4( const double x , void *userdata );
/* 1d cquad called by IntegrandCC_ndiag */
double IntegrandCC_ndiag  ( const double x,  void *userdata );
double Integrand1d_cquadn1( const double x , void *userdata );
double Integrand1d_cquadn2( const double x , void *userdata );
double Integrand1d_cquadn3( const double x , void *userdata );
double Integrand1d_cquadn4( const double x , void *userdata );
double Integrand1d_cquadn5( const double x , void *userdata );
double Integrand1d_cquadn6( const double x , void *userdata );
}


#endif /* SCATTERING_H */

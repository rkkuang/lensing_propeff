static char help[] = "Basic Lensing operator and Scattering operator.\n";
/*
 
 Example code of lensing with scattering effect
 
 This file can be compiled with "make scat" using the example PETSc makefile at 
 https://github.com/petsc/petsc/blob/main/share/petsc/Makefile.user

 There are three output *.dat files corresponding to the source, lensed image, and lensed image with scattering effect, which are in ASCII format, and can be visualized with pyplot.py   

 */

#include <petscmat.h>
#include <time.h>

// functions
void exp_brightness(PetscReal scale_length, PetscReal central_bright, PetscReal ax_ratio, PetscReal PA, PetscReal physical_centerx, PetscReal physical_centery, PetscInt size, PetscReal *xs, PetscReal *ys, PetscReal *Iprof);
void pixel_physical_position(PetscReal resolution, PetscInt xpixNum, PetscInt ypixNum, PetscReal img_centerx, PetscReal img_centery, PetscInt size, PetscInt *II, PetscInt *JJ, PetscReal corrx, PetscReal corry, PetscReal *physx, PetscReal *physy);
void SIE_deflect_ang(PetscInt size, PetscReal *xs, PetscReal *ys, PetscReal re, PetscReal Lxc, PetscReal Lyc, PetscReal PA, PetscReal f, PetscReal *dx, PetscReal *dy);
void sub_setMtx(PetscReal px, PetscReal py, PetscReal resolution, PetscReal img_centerx, PetscReal img_centery, PetscInt L, PetscInt j, PetscReal const1, PetscReal const2, Mat Mtx);
void getLmtx(PetscInt size, PetscReal *physx_src, PetscReal *physy_src, PetscReal *ext, PetscReal resolution, PetscReal img_centerx, PetscReal img_centery, PetscInt *mnkl, PetscReal corrx, PetscReal corry, Mat Lmtx);
void getSmtx(PetscInt size, PetscReal *physx_src, PetscReal *physy_src, PetscReal *ext, PetscReal resolution, PetscReal img_centerx, PetscReal img_centery, PetscInt *mnkl, PetscReal corrx, PetscReal corry, Mat Smtx, PetscInt *scat_idx, PetscInt scat_idx_cnt);
void genxyint(PetscInt xpixNum, PetscInt ypixNum, PetscInt *II, PetscInt *JJ);


// source surface brightness profile
void exp_brightness(PetscReal scale_length, PetscReal central_bright, PetscReal ax_ratio, PetscReal PA, PetscReal physical_centerx, PetscReal physical_centery, PetscInt size, PetscReal *xs, PetscReal *ys, PetscReal *Iprof)
{
  PetscReal costheta, sintheta, xs_Xc, ar2, ys_Yc, u, v, dis;
  PetscInt i;
  costheta = cos(PA);
  sintheta = sin(PA);
  ar2 = ax_ratio * ax_ratio;

  for (i = 0; i < size; i++) {
    xs_Xc = xs[i] - physical_centerx;
    ys_Yc = ys[i] - physical_centery;
    u = xs_Xc * costheta + ys_Yc * sintheta;
    v = -xs_Xc * sintheta + ys_Yc * costheta;

    dis = sqrt( u * u / ar2 + v * v );

    Iprof[i] = central_bright * exp(-dis / scale_length);
  }
}

// given a pixel coordinate, return the physical coordinate in e.g. arcsec
void pixel_physical_position(PetscReal resolution, PetscInt xpixNum, PetscInt ypixNum, PetscReal img_centerx, PetscReal img_centery, PetscInt size, PetscInt *II, PetscInt *JJ, PetscReal corrx, PetscReal corry, PetscReal *physx, PetscReal *physy) {
// '''
// resolution: arcsec per pixel
// xpixNum, ypixNum: how many pixels in each side
// img_center: the physical position (in arcsec) of the image center

// return: physical positions of image pixel centers
// '''
  PetscInt i;
  for (i = 0; i < size; i++) {
    physx[i] = img_centerx + (JJ[i] + corrx - xpixNum * 0.5) * resolution;
    physy[i] = img_centery + (ypixNum * 0.5 - II[i] + corry) * resolution;
  }
}


// calculate the deflection angle of a SIE lens
void SIE_deflect_ang(PetscInt size, PetscReal *xs, PetscReal *ys, PetscReal re, PetscReal Lxc, PetscReal Lyc, PetscReal PA, PetscReal f, PetscReal *dx, PetscReal *dy)
{
  PetscReal cosPA, sinPA, fprime, f_sqrt, coe_farsinh , coe_fcosphi, xs_Xc, ys_Yc, u, v, phis, cosphis, sinphis, newx, newy, du, dv;
  PetscInt i;

  // https://petsc.org/release/include/petscmath.h.html
  cosPA = cos(PA);
  sinPA = sin(PA);
  fprime = sqrt(1 - f * f);
  f_sqrt = sqrt(f);
  coe_farsinh = f_sqrt / fprime;
  coe_fcosphi = fprime / f;

  for (i = 0; i < size; i++) {
    xs_Xc = xs[i] - Lxc;
    ys_Yc = ys[i] - Lyc;
    u = xs_Xc * cosPA + ys_Yc * sinPA;
    v = -xs_Xc * sinPA + ys_Yc * cosPA;

    phis = atan2(v, u);
    cosphis = cos(phis);
    sinphis = sin(phis);

    du = -coe_farsinh * asinh(coe_fcosphi * cosphis) * re;
    dv = -coe_farsinh * asin(fprime * sinphis) * re;
    u += du;
    v += dv;

    newx = u * cosPA - v * sinPA + Lxc;
    newy = u * sinPA + v * cosPA + Lyc;

    dx[i] = xs[i] - newx;
    dy[i] = ys[i] - newy;

  }
}


void sub_setMtx(PetscReal px, PetscReal py, PetscReal resolution, PetscReal img_centerx, PetscReal img_centery, PetscInt L, PetscInt j, PetscReal const1, PetscReal const2, Mat Mtx) {
PetscReal px_pix, py_pix, corner1x, corner1y, corner2x, corner2y, xweight, yweight;
PetscInt xcorrect, ycorrect, jj, I1, J1, I2, J2, I3, J3, I4, J4;   
PetscScalar    v;   
      //# determine which 4 pixels in source plane enclose (px, py), and the corresponding weights, then place them in the j-th raw of the Lmtx array
      px -= img_centerx;
      py -= img_centery;
      px_pix = (px / resolution);
      py_pix = (py / resolution);

      xcorrect = (px_pix > 0) ? -1 : 1;
      ycorrect = (py_pix > 0) ? -1 : 1;

      corner1x = (int)(round(px_pix)) + 0.5 * xcorrect;
      corner1y = (int)(round(py_pix)) + 0.5 * ycorrect;

      yweight = fabs(py_pix - corner1y);
      xweight = fabs(px_pix - corner1x);

      corner2x = corner1x - xcorrect;
      corner2y = corner1y - ycorrect;

      J1 = (int)(-(corner1y + const1)); //# numpy 数组 里的行
      I1 = (int)(corner1x + const2); //# numpy 数组里的列
      jj = I1 * L + J1; v = (1 - xweight) * (1 - yweight);
      MatSetValues(Mtx, 1, &j, 1, &jj, &v, INSERT_VALUES);

      J2 = (int)(-(corner2y + const1));
      I2 = (int)(corner1x + const2);
      jj = I2 * L + J2; v = (1 - xweight) * yweight;
      MatSetValues(Mtx, 1, &j, 1, &jj, &v, INSERT_VALUES);

      J3 = (int)(-(corner2y + const1));
      I3 = (int)(corner2x + const2);
      jj = I3 * L + J3; v = xweight * yweight;
      MatSetValues(Mtx, 1, &j, 1, &jj, &v, INSERT_VALUES);

      J4 = (int)(-(corner1y + const1));
      I4 = (int)(corner2x + const2);
      jj = I4 * L + J4; v = xweight * (1 - yweight);
      MatSetValues(Mtx, 1, &j, 1, &jj, &v, INSERT_VALUES);

}

// construct the Lensing operator
void getLmtx(PetscInt size, PetscReal *physx_src, PetscReal *physy_src, PetscReal *ext, PetscReal resolution, PetscReal img_centerx, PetscReal img_centery, PetscInt *mnkl, PetscReal corrx, PetscReal corry, Mat Lmtx) {
// # we should obtain the Lensing operator:
  PetscInt M, N, K, L, xpixNum, ypixNum, i, row_start, row_end;//, xcorrect, ycorrect, I1, J1, I2, J2, I3, J3, I4, J4;
  // PetscBool cnd1, cnd2, cnd3, cnd4, cnd;
  PetscReal const1, const2, px, py;// px_pix, py_pix, corner1x, corner1y, corner2x, corner2y, xweight, yweight;

  M = mnkl[0]; N = mnkl[1]; K = mnkl[2]; L = mnkl[3];
  xpixNum = K; ypixNum = L;
  //Lmtx = np.zeros((M*N, K*L)) # 3600 x 900

  const1 = -corry - ypixNum * 0.5;
  const2 = -corrx + xpixNum * 0.5;
  MatZeroEntries(Lmtx);
  MatGetOwnershipRange(Lmtx, &row_start, &row_end);

  for (i = row_start; i < row_end; i++){

    if ( (physx_src[i] > (ext[0] + resolution)) && (physx_src[i] < (ext[1] - resolution)) && (physy_src[i] > (ext[2] + resolution)) && (physy_src[i] < (ext[3] - resolution))) {
      px = physx_src[i]; py = physy_src[i];
      sub_setMtx(px, py, resolution, img_centerx, img_centery, L, i, const1, const2, Lmtx);
    }

  }

  // for (i = 0; i < size; i++) {
  //   if ( (physx_src[i] > (ext[0] + resolution)) && (physx_src[i] < (ext[1] - resolution)) && (physy_src[i] > (ext[2] + resolution)) && (physy_src[i] < (ext[3] - resolution))) {
  //     px = physx_src[i]; py = physy_src[i];
  //     sub_setMtx(px, py, resolution, img_centerx, img_centery, L, i, const1, const2, Lmtx);
  //   }
  // }

  MatAssemblyBegin(Lmtx, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Lmtx, MAT_FINAL_ASSEMBLY);
}


void getSmtx(PetscInt size, PetscReal *physx_src, PetscReal *physy_src, PetscReal *ext, PetscReal resolution, PetscReal img_centerx, PetscReal img_centery, PetscInt *mnkl, PetscReal corrx, PetscReal corry, Mat Smtx, PetscInt *scat_idx, PetscInt scat_idx_cnt) {
  PetscInt M, N, K, L, xpixNum, ypixNum, i, ii, row_start, row_end;//, j, ii, jj, xcorrect, ycorrect, I1, J1, I2, J2, I3, J3, I4, J4;
  PetscReal const1, const2, px, py;//, px_pix, py_pix, corner1x, corner1y, corner2x, corner2y, xweight, yweight;
  PetscScalar    v;
  M = mnkl[0]; N = mnkl[1]; K = mnkl[2]; L = mnkl[3];
  xpixNum = K; ypixNum = L;
  const1 = -corry - ypixNum * 0.5;
  const2 = -corrx + xpixNum * 0.5;

  MatZeroEntries(Smtx);


  // for (ii = 0; ii < size; ii++) {
  //   v = 1.0;
  //   MatSetValues(Smtx, 1, &ii, 1, &ii, &v, INSERT_VALUES);
  // }

  MatGetOwnershipRange(Smtx, &row_start, &row_end);
  v = 1.0;
  for (ii = row_start; ii < row_end; ii++){
    MatSetValues(Smtx, 1, &ii, 1, &ii, &v, INSERT_VALUES);
  }


  for (ii = 0; ii < scat_idx_cnt; ii++) {
    i = scat_idx[ii];
    if ( (physx_src[i] > (ext[0] + resolution)) && (physx_src[i] < (ext[1] - resolution)) && (physy_src[i] > (ext[2] + resolution)) && (physy_src[i] < (ext[3] - resolution))) {
      v = 0;
      MatSetValues(Smtx, 1, &i, 1, &i, &v, INSERT_VALUES);

      px = physx_src[i]; py = physy_src[i];
      sub_setMtx(px, py, resolution, img_centerx, img_centery, L, i, const1, const2, Smtx);
    }
  }
  MatAssemblyBegin(Smtx, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Smtx, MAT_FINAL_ASSEMBLY);
}

void genxyint(PetscInt xpixNum, PetscInt ypixNum, PetscInt *II, PetscInt *JJ) {
  PetscInt i, j, k;
  for (i = 0; i < xpixNum; i++) {
    k = i * ypixNum;
    for (j = 0; j < ypixNum; j++) {
      II[k + j] = j;
      JJ[k + j] = i;
    }
  }
}




// // https://stackoverflow.com/questions/7034930/how-to-generate-gaussian-pseudo-random-numbers-in-c-for-a-given-mean-and-varianc
// double drand()   /* uniform distribution, (0..1] */
// {
//   return (rand() + 1.0) / (RAND_MAX + 1.0);
// }

// double random_normal()
// /* normal distribution, centered on 0, std dev 1 */
// {
//   return sqrt(-2 * log(drand())) * cos(2 * M_PI * drand());
// }

/*
 * Normal random numbers generator - Marsaglia algorithm.
 */
double* generate(int n)
{

  int i;
  int m = n + n % 2;
  double* values = (double*)calloc(m, sizeof(double));
  if ( values )
  {
    for ( i = 0; i < m; i += 2 )
    {
      double x, y, rsq, f;
      do {
        x = 2.0 * rand() / (double)RAND_MAX - 1.0;
        y = 2.0 * rand() / (double)RAND_MAX - 1.0;
        rsq = x * x + y * y;
      } while ( rsq >= 1. || rsq == 0. );
      f = sqrt( -2.0 * log(rsq) / rsq );
      values[i]   = x * f;
      values[i + 1] = y * f;
    }
  }
  return values;
}


int main(int argc, char **argv)
{
  PetscReal corrx, corry, resolution, resol_imgplane, xyregion, img_centerx, img_centery;// img_plane_centerx, img_plane_centery;
  PetscInt xpixNum, ypixNum,  M, N, Nsrc, Nimg, i, scat_times, nlimits = 4;
  PetscReal physical_center1x, physical_center1y, scale_length1, central_bright1, ax_ratio1, PA1;
  PetscReal physical_center2x, physical_center2y, scale_length2, central_bright2, ax_ratio2, PA2;




  corrx = 0.5; corry = -0.5;

  // # source plane grid number and resolution
  // xpixNum = 60;
  // ypixNum = 60; 
  
  

  // # image plane pixel number and resolution
  
  PetscReal sieThetaE, siecx, siecy, siePA, siear;

  // PetscInitialize(&argc, &argv, (char*)0, help);

  // PetscInitialize(int *argc, char ***args,const char *file,const char *help); /* C */

  PetscReal scat_reg[4];// = { -0.2, 1.2, -1.25, 0.2};

  PetscReal scat_angle;
  // scat_angle = 0.05;

  PetscInitialize(&argc, &argv, "pOption.yaml", help);

  PetscOptionsGetInt(NULL,NULL,"-pixnumSrc",&ypixNum, NULL);
  PetscOptionsGetInt(NULL,NULL,"-pixnumImg",&M, NULL);

  PetscOptionsGetReal(NULL,NULL,"-regSrc",&resolution, NULL);
  PetscOptionsGetReal(NULL,NULL,"-regImg",&resol_imgplane, NULL);
  PetscOptionsGetReal(NULL,NULL,"-imgcx",&img_centerx, NULL);
  PetscOptionsGetReal(NULL,NULL,"-imgcy",&img_centery, NULL);

  PetscOptionsGetReal(NULL,NULL,"-src1cx",&physical_center1x, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src1cy",&physical_center1y, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src1r",&scale_length1, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src1I",&central_bright1, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src1ar",&ax_ratio1, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src1PA",&PA1, NULL);

  PetscOptionsGetReal(NULL,NULL,"-src2cx",&physical_center2x, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src2cy",&physical_center2y, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src2r",&scale_length2, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src2I",&central_bright2, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src2ar",&ax_ratio2, NULL);
  PetscOptionsGetReal(NULL,NULL,"-src2PA",&PA2, NULL);

  PetscOptionsGetReal(NULL,NULL,"-sieThetaE",&sieThetaE, NULL);
  PetscOptionsGetReal(NULL,NULL,"-siecx",&siecx, NULL);
  PetscOptionsGetReal(NULL,NULL,"-siecy",&siecy, NULL);
  PetscOptionsGetReal(NULL,NULL,"-siePA",&siePA, NULL);
  PetscOptionsGetReal(NULL,NULL,"-siear",&siear, NULL);
  
  PetscOptionsGetReal(NULL,NULL,"-scat_angle",&scat_angle, NULL);
  PetscOptionsGetRealArray(NULL,NULL,"-scat_reg", scat_reg, &nlimits,  NULL);
  PetscOptionsGetInt(NULL,NULL,"-rand_trails",&scat_times, NULL);




  // // // ### build the source vector
  // physical_center1x = -0.05; physical_center1y = 0.05;
  // scale_length1 = 0.1;
  // central_bright1 = 100;
  // ax_ratio1 = 0.64;
  PA1 = PA1 / 180.0 * M_PI;

  // physical_center2x = -0.4; physical_center2y = 0.25;
  // scale_length2 = 0.1;
  // central_bright2 = 50;
  // ax_ratio2 = 1;
  PA2 = PA2 / 180.0 * M_PI;


  // PetscReal re, PetscReal Lxc, PetscReal Lyc, PetscReal PA, PetscReal f
  // 0.9, 0, 0, 45.0 / 180 * M_PI, 0.8




  // ./scat -pixnumSrc 60 -pixnumImg 100 -regSrc 1.5 -regImg 4 -imgcx -0.17 -imgcy 0.2

  
  xpixNum = ypixNum;
  resolution = resolution / ypixNum;
  Nsrc = xpixNum * ypixNum;
  N = M; Nimg = M * N;

  xyregion = resol_imgplane;
  resol_imgplane = resol_imgplane / M; // #xyregion/M
  // img_centerx = -0.17; img_centery = 0.2;

  //xyregion = M * resol_imgplane;




  srand(time(NULL));

  PetscReal ext[4] = { -resolution*(xpixNum / 2) + img_centerx, resolution*(xpixNum / 2) + img_centerx, -resolution*(ypixNum / 2) + img_centery, resolution*(ypixNum / 2) + img_centery};

  PetscInt II[Nsrc], JJ[Nsrc], II2[Nimg], JJ2[Nimg];
  PetscReal physx[Nsrc], physy[Nsrc], physx_lens[Nimg], physy_lens[Nimg], Iprof1[Nsrc], Iprof2[Nsrc], dxSIE[Nimg], dySIE[Nimg], physx_src[Nimg], physy_src[Nimg];// physx_src2[Nimg], physy_src2[Nimg];

  genxyint(xpixNum, ypixNum, II, JJ);





  pixel_physical_position(resolution, xpixNum, ypixNum, img_centerx, img_centery, Nsrc, II, JJ, corrx, corry, physx, physy);

  exp_brightness(scale_length1, central_bright1, ax_ratio1, PA1, physical_center1x, physical_center1y, Nsrc, physx, physy, Iprof1);
  exp_brightness(scale_length2, central_bright2, ax_ratio2, PA2, physical_center2x, physical_center2y, Nsrc, physx, physy, Iprof2);

  Vec src1d;
  PetscScalar tmp_scalar;
  VecCreate(PETSC_COMM_WORLD, &src1d); // In most cases users can employ the communicator PETSC_COMM_WORLD to indicate all processes in a given run and PETSC_COMM_SELF to indicate a single process.
  VecSetSizes(src1d, PETSC_DECIDE, Nsrc);
  VecSetFromOptions(src1d);

  // fprintf(stderr, "can you 245 \n");

  for (int i = 0; i < Nsrc; i++) {
    tmp_scalar = Iprof1[i] + Iprof2[i];
    VecSetValue(src1d, i, tmp_scalar, INSERT_VALUES);
  }
  VecAssemblyBegin(src1d);
  VecAssemblyEnd(src1d);

  PetscScalar tmp_scalar2[5];
  PetscInt tmp_scalar2_idx[5] = {0, 1, 2, 3, 4};

  VecGetValues(src1d, 5, tmp_scalar2_idx, tmp_scalar2);

  PetscReal extimg[4] = { -xyregion / 2 + img_centerx, xyregion / 2 + img_centerx, -xyregion / 2 + img_centery, xyregion / 2 + img_centery};

  // # lensing operator: Lmtx, size MN raws and KL columns, i.e., 3600 x 900
  // # first, each pixel j = m + (n-1)*M in image plane is cast back on the source plane to a position yj, using the lens equation
  // # obtain the physical position of all pixels in the lens plane

  genxyint(M, N, II2, JJ2);

  pixel_physical_position(resol_imgplane, M, N, img_centerx, img_centery, Nimg, II2, JJ2, corrx, corry, physx_lens, physy_lens);

  // # using the lens equation, cast images to the source plane
  SIE_deflect_ang(Nimg, physx_lens, physy_lens, sieThetaE, siecx, siecy, siePA / 180 * M_PI, siear, dxSIE, dySIE);
  //sieThetaE, siecx, siecy, siePA, siear;

  // PetscReal scat_angle, scatx[Nimg], scaty[Nimg], scatx_reg[Nimg] = 0, scaty_reg[Nimg] = 0, physx_lens_[Nimg], physy_lens_[Nimg];
  PetscReal scatx_reg, scaty_reg, physx_lens_[Nimg], physy_lens_[Nimg];



  PetscInt scat_idx[Nimg], scat_idx_cnt = 0;

  for (int i = 0; i < Nimg; i++) {
    physx_src[i] = physx_lens[i] - dxSIE[i];
    physy_src[i] = physy_lens[i] - dySIE[i];
  }

  Mat Lmtx, Smtx, Smtx_sum; //Lmtx2
  MatCreate(PETSC_COMM_WORLD, &Lmtx); // PETSC_COMM_WORLD -- multi process
  MatSetSizes(Lmtx, PETSC_DECIDE, PETSC_DECIDE, Nimg, Nsrc);
  // MatSetType(Lmtx, MATAIJ);
  MatSetFromOptions(Lmtx);
  // MatSetUp(Lmtx);


  // MatMPIAIJSetPreallocation(Lmtx, 2, NULL, 2, NULL);
  // Preallocates memory for a sparse parallel matrix in AIJ format (the default parallel PETSc format).
  MatSeqAIJSetPreallocation(Lmtx, 4, NULL);


  // MatGetOwnershipRange(Lmtx, &row_start, &row_end);
  // Returns the range of matrix rows owned by this processor, assuming that the matrix is laid out with the first n1 rows on the first processor, the next n2 rows on the second, etc. For certain parallel layouts this range may not be well defined.

  // m - the global index of the first local row
  // n - one more than the global index of the last local row



  PetscInt mnkl[4] = {M, N, xpixNum, ypixNum};
  PetscInt mnmn[4] = {M, N, M, N};

  getLmtx(Nimg, physx_src, physy_src, ext, resolution, img_centerx, img_centery, mnkl, corrx, corry, Lmtx);

  Vec img1d, img1d_scat, img1d_scat_sum; //img1d2
  VecCreate(PETSC_COMM_WORLD, &img1d);
  VecSetSizes(img1d, PETSC_DECIDE, Nimg);
  VecSetFromOptions(img1d);

  VecCreate(PETSC_COMM_WORLD, &img1d_scat);
  VecSetSizes(img1d_scat, PETSC_DECIDE, Nimg);
  VecSetFromOptions(img1d_scat);

  VecCreate(PETSC_COMM_WORLD, &img1d_scat_sum);
  VecSetSizes(img1d_scat_sum, PETSC_DECIDE, Nimg);
  VecSetFromOptions(img1d_scat_sum);

  MatMult(Lmtx, src1d, img1d);

  MatCreate(PETSC_COMM_WORLD, &Smtx);
  MatSetSizes(Smtx, PETSC_DECIDE, PETSC_DECIDE, Nimg, Nimg);
  MatSetType(Smtx, MATAIJ);
  MatSetFromOptions(Smtx);
  // MatSeqAIJSetPreallocation(Smtx, 4, NULL);
  MatSetUp(Smtx);

  MatCreate(PETSC_COMM_WORLD, &Smtx_sum);
  MatSetSizes(Smtx_sum, PETSC_DECIDE, PETSC_DECIDE, Nimg, Nimg);
  MatSetType(Smtx_sum, MATAIJ);
  MatSetFromOptions(Smtx_sum);
  MatSetUp(Smtx_sum);


  PetscInt k, ii;
  PetscScalar vecsum;
  PetscReal matnorm;

  


  // Computes Y = a*X + Y.
  // MatAXPY(Mat Y,PetscScalar a,Mat X,MatStructure str) //https://petsc.org/release/docs/manualpages/Mat/MatAXPY.html

  for (i = 0; i < Nimg; i++) {
    if (physx_lens[i] > scat_reg[0] && physx_lens[i] < scat_reg[1] && physy_lens[i] > scat_reg[2] && physy_lens[i] < scat_reg[3]) {
      scat_idx[scat_idx_cnt] = i;
      scat_idx_cnt += 1;
    }
    physx_lens_[i] = physx_lens[i];
    physy_lens_[i] = physy_lens[i];
  }

  VecSet(img1d_scat_sum, 0);
  MatZeroEntries(Smtx_sum);
  MatAssemblyBegin(Smtx_sum, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Smtx_sum, MAT_FINAL_ASSEMBLY);

  double *randxy;
  for (k = 0; k < scat_times; k++) {

    for (ii = 0; ii < scat_idx_cnt; ii++) {
      i = scat_idx[ii];
      randxy = generate(2);

      scatx_reg = 0 + scat_angle * randxy[0];//random_normal();
      scaty_reg = 0 + scat_angle * randxy[1];//random_normal();

      physx_lens_[i] = physx_lens[i] + scatx_reg;
      physy_lens_[i] = physy_lens[i] + scaty_reg;

    }

    
    getSmtx(Nimg, physx_lens_, physy_lens_, extimg, resol_imgplane, img_centerx, img_centery, mnmn, corrx, corry, Smtx, scat_idx, scat_idx_cnt);

    // MatMult(Smtx, img1d, img1d_scat); // scattered image

    // VecAXPY(img1d_scat_sum, 1, img1d_scat);
    
    MatAXPY(Smtx_sum, 1.0, Smtx, DIFFERENT_NONZERO_PATTERN);

  }

  // VecScale(img1d_scat_sum, 1.0 / scat_times);
  


  MatScale(Smtx_sum, 1.0 / scat_times);



  MatMult(Smtx_sum, img1d, img1d_scat_sum); // scattered image
  

  VecSum(img1d, &vecsum);
  printf("sum w/o scattering %f \n", vecsum);
  VecSum(img1d_scat_sum, &vecsum);
  printf("sum w/ scattering %f \n", vecsum);

  VecSum(src1d, &vecsum);
  printf("sum source %f \n", vecsum);

  MatNorm(Lmtx, NORM_FROBENIUS, &matnorm); // NORM_1, NORM_FROBENIUS
  printf("Lmtx L2 norm = %f \n", matnorm);

  MatNorm(Smtx_sum, NORM_FROBENIUS, &matnorm);
  printf("Smtx L2 norm = %f \n", matnorm);


  // save out src1d, img1d and img1d_scat

  PetscViewer viewer;

  /* open ASCII file for writing */
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);

  PetscViewerFileSetName(viewer, "./src1d.dat");
  VecView(src1d, viewer);

  PetscViewerFileSetName(viewer, "./img1d.dat");
  VecView(img1d, viewer);

  PetscViewerFileSetName(viewer, "./img1d_scat.dat");
  VecView(img1d_scat_sum, viewer);

  PetscViewerFileSetName(viewer, "./Smtx.dat");
  MatView(Smtx_sum, viewer);

  PetscViewerFileSetName(viewer, "./Lmtx.dat");
  MatView(Lmtx, viewer);



  PetscViewerDestroy(&viewer);


  VecDestroy(&src1d);
  VecDestroy(&img1d);
  VecDestroy(&img1d_scat);
  VecDestroy(&img1d_scat_sum);
  MatDestroy(&Lmtx);
  MatDestroy(&Smtx);
  MatDestroy(&Smtx_sum);


  PetscFinalize();
  return 0;
}
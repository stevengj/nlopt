/*
Copyright (C) 2018 Sovrasov V. - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license.
 * You should have received a copy of the MIT license with
 * this file. If not visit https://opensource.org/licenses/MIT
*/
#include "evolvent.hpp"
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace ags;

void mapd(double x, int m, double* y, int n, int key = 1);    /* map x to y         */
void invmad(int, double [], int, int *, double [], int, int);  /* map y to x         */
void xyd(double *xx, int m, double y[], int n);        /* get preimage       */

Evolvent::Evolvent()
{
  mIsInitialized = false;
}

Evolvent::~Evolvent()
{
}

Evolvent::Evolvent(int dimension, int tightness, const double* lb, const double* ub, MapType type)
{
  assert(tightness > 2);

  mDimension = dimension;
  mTightness = tightness;
  mMapType = type;

  p2.resize(mDimension);
  mShiftScalars.resize(mDimension);
  mRho.resize(mDimension);
  for (int i = 0; i < mDimension; i++)
  {
    mRho[i] = ub[i] - lb[i];
    mShiftScalars[i] = 0.5*(lb[i] + ub[i]);
  }

  switch (mMapType)
  {
  case Simple:
    mMapKey = 1;
    break;
  case Linear:
    mMapKey = 2;
    break;
  case Noninjective:
    mMapKey = 3;
    break;
  }

  mIsInitialized = true;
}

void Evolvent::TransformToStandardCube(const double *y, double *z)
{
  for (int i = 0; i < mDimension; i++)
    z[i] = (y[i] - mShiftScalars[i]) / mRho[i];
}

void Evolvent::TransformToSearchDomain(const double *y, double *z)
{
  for (int i = 0; i < mDimension; i++)
    z[i] = mRho[i]*y[i] + mShiftScalars[i];
}

void Evolvent::GetImage(double x, double y[])
{
  if(mDimension != 1)
    mapd(x, mTightness, y, mDimension, mMapKey);
  else
    y[0] = x - 0.5;

  TransformToSearchDomain(y, y);
}

int Evolvent::GetAllPreimages(const double *p, double xp[])
{
  int preimNumber = 1;
  TransformToStandardCube(p, p2.data());
  if(mMapType == Noninjective)
    invmad(mTightness, xp, noninjectiveMaxPreimages, &preimNumber, p2.data(), mDimension, 4);
  else
    xyd(xp, mTightness, p2.data(), mDimension);

  return preimNumber;
}

void xyd(double *xx, int m, double y[], int n)
{
  /* calculate preimage x  for nearest level  m center to y */
  /* (x - left boundary point of level m interval)          */
  int n1, nexp, l, iu[10], iv[10];

  double x, r1;
  double r;
  int iw[11];
  int i, j, it, is;
  void numbr(int *, const int, const int, int&, int*, int*);

  n1 = n - 1;
  for (nexp = 1, i = 0; i<n; i++) {
    nexp *= 2; iw[i] = 1;
  }
  r = 0.5;
  r1 = 1.0;
  x = 0.0;
  it = 0;
  for (j = 0; j<m; j++) {
    r *= 0.5;
    for (i = 0; i<n; i++) {
      iu[i] = (y[i]<0) ? -1 : 1;
      y[i] -= r*iu[i];
      iu[i] *= iw[i];
    }
    i = iu[0];
    iu[0] = iu[it];
    iu[it] = i;
    numbr(&is, n1, nexp, l, iu, iv);
    i = iv[0];
    iv[0] = iv[it];
    iv[it] = i;
    for (i = 0; i<n; i++)
      iw[i] = -iw[i] * iv[i];
    if (l == 0) l = it;
    else if (l == it) l = 0;
    it = l;
    r1 = r1 / nexp;
    x += r1*is;
  }
  *xx = x;
}
void numbr(int *iss, const int n1, const int nexp, int& l, int* iu, int* iv)
{
  /* calculate s(u)=is,l(u)=l,v(u)=iv by u=iu */

  int i, n, is, iff, k1, k2, l1;

  n = n1 + 1;
  iff = nexp;
  is = 0;
  k1 = -1;
  for (i = 0; i<n; i++) {
    iff = iff / 2;
    k2 = -k1*iu[i];
    iv[i] = iu[i];
    k1 = k2;
    if (k2<0) l1 = i;
    else { is += iff; l = i; }
  }
  if (is == 0) l = n1;
  else {
    iv[n1] = -iv[n1];
    if (is == (nexp - 1)) l = n1;
    else if (l1 == n1) iv[l] = -iv[l];
    else l = l1;
  }
  *iss = is;
}

void invmad(int m, double xp[], int kp,
  int *kxx, double p[], int n, int incr)
{
  /* calculate kx preimage p node */
  /*   node type mapping m level  */

  double mne, d1, dd, x, dr, del;// , convers;
  double r, d, u[10], y[10];
  int i, k, kx, nexp;
  void xyd(double *, int, double *, int);

  kx = 0;
  kp--;
  for (nexp = 1, i = 0; i<n; i++) {
    nexp *= 2; u[i] = -1.0;
  }
  dr = nexp;
  for (mne = 1, r = 0.5, i = 0; i<m; i++) {
    mne *= dr; r *= 0.5;
  }
  dr = mne / nexp;

  dr = dr - fmod(dr, 1.0);
  //dr = (dr>0) ? floor(dr) : ceil(dr);

  del = 1. / (mne - dr);
  d1 = del*(incr + 0.5);
  for (kx = -1; kx<kp;) {
    for (i = 0; i<n; i++) {       /* label 2 */
      d = p[i];
      y[i] = d - r*u[i];
    }
    for (i = 0; (i<n) && (fabs(y[i]) < 0.5); i++);
    if (i >= n) {
      xyd(&x, m, y, n);
      dr = x*mne;
      dd = dr - fmod(dr, 1.0);
      //dd = (dr>0) ? floor(dr) : ceil(dr);
      dr = dd / nexp;
      dd = dd - dr + fmod(dr, 1.0);
      //dd = dd - ((dr>0) ? floor(dr) : ceil(dr));
      x = dd*del;
      if (kx>kp) break;
      k = kx++;                     /* label 9 */
      if (kx == 0) xp[0] = x;
      else {
        while (k >= 0) {
          dr = fabs(x - xp[k]);     /* label 11 */
          if (dr <= d1) {
            for (kx--; k<kx; k++, xp[k] = xp[k + 1]);
            goto m6;
          }
          else
            if (x <= xp[k]) {
              xp[k + 1] = xp[k]; k--;
            }
            else break;
        }
        xp[k + 1] = x;
      }
    }
  m6: for (i = n - 1; (i >= 0) && (u[i] = (u[i] <= 0.0) ? 1 : -1)<0; i--);
  if (i<0)break;
  }
  *kxx = ++kx;

}

void mapd(double x, int m, double* y, int n, int key)
{
  /* mapping y(x) : 1 - center, 2 - line, 3 - node */
  // use key = 1

  int n1, nexp, l, iq, iu[10], iv[10];
  double d, mne, dd, dr;//,tmp;
  double p, r;
  int iw[11];
  int it, is, i, j, k;
  void node(int is, int n1, int nexp, int& l, int& iq, int iu[], int iv[]);

  p = 0.0;
  n1 = n - 1;
  for (nexp = 1, i = 0; i<n; nexp *= 2, i++);   /* nexp=2**n */
  d = x;
  r = 0.5;
  it = 0;
  dr = nexp;
  for (mne = 1, i = 0; i<m; mne *= dr, i++);    /* mne=dr**m  */
  for (i = 0; i<n; i++) {
    iw[i] = 1; y[i] = 0.0;
  }

  if (key == 2) {
    d = d*(1.0 - 1.0 / mne); k = 0;
  }
  else
    if (key > 2) {
      dr = mne / nexp;
      dr = dr - fmod(dr, 1.0);
      //dr=(dr>0)?floor(dr):ceil(dr);
      dd = mne - dr;
      dr = d*dd;
      dd = dr - fmod(dr, 1.0);
      //dd=(dr>0)?floor(dr):ceil(dr);
      dr = dd + (dd - 1) / (nexp - 1);
      dd = dr - fmod(dr, 1.0);
      //dd=(dr>0)?floor(dr):ceil(dr);
      d = dd*(1. / (mne - 1.0));
    }

  for (j = 0; j<m; j++) {
    iq = 0;
    if (x == 1.0) {
      is = nexp - 1; d = 0.0;
    }
    else {
      d = d*nexp;
      is = (int)d;
      d = d - is;
    }
    i = is;
    node(i, n1, nexp, l, iq, iu, iv);
    i = iu[0];
    iu[0] = iu[it];
    iu[it] = i;
    i = iv[0];
    iv[0] = iv[it];
    iv[it] = i;
    if (l == 0)
      l = it;
    else if (l == it) l = 0;
    if ((iq>0) || ((iq == 0) && (is == 0)))  k = l;
    else if (iq<0) k = (it == n1) ? 0 : n1;
    r = r*0.5;
    it = l;
    for (i = 0; i<n; i++) {
      iu[i] = iu[i] * iw[i];
      iw[i] = -iv[i] * iw[i];
      p = r*iu[i];
      p = p + y[i];
      y[i] = p;
    }
  }
  if (key == 2) {
    if (is == (nexp - 1)) i = -1;
    else i = 1;
    p = 2 * i*iu[k] * r*d;
    p = y[k] - p;
    y[k] = p;
  }
  else if (key == 3) {
    for (i = 0; i<n; i++) {
      p = r*iu[i];
      p = p + y[i];
      y[i] = p;
    }
  }
}

void node(int is, int n1, int nexp, int& l, int& iq, int iu[], int iv[])
{
  /* calculate iu=u[s], iv=v[s], l=l[s] by is=s */

  int n, i, j, k1, k2, iff;

  n = n1 + 1;
  if (is == 0) {
    l = n1;
    for (i = 0; i<n; i++) {
      iu[i] = -1; iv[i] = -1;
    }
  }
  else if (is == (nexp - 1)) {
    l = n1;
    iu[0] = 1;
    iv[0] = 1;
    for (i = 1; i<n; i++) {
      iu[i] = -1; iv[i] = -1;
    }
    iv[n1] = 1;
  }
  else {
    iff = nexp;
    k1 = -1;
    for (i = 0; i<n; i++) {
      iff = iff / 2;
      if (is >= iff) {
        if ((is == iff) && (is != 1)) { l = i; iq = -1; }
        is = is - iff;
        k2 = 1;
      }
      else {
        k2 = -1;
        if ((is == (iff - 1)) && (is != 0)) { l = i; iq = 1; }
      }
      j = -k1*k2;
      iv[i] = j;
      iu[i] = j;
      k1 = k2;
    }
    iv[l] = iv[l] * iq;
    iv[n1] = -iv[n1];
  }
}

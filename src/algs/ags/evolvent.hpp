/*
Copyright (C) 2018 Sovrasov V. - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license.
 * You should have received a copy of the MIT license with
 * this file. If not visit https://opensource.org/licenses/MIT
*/
#pragma once

#include <vector>

namespace ags
{

enum MapType {
  Simple = 1, Linear = 2, Noninjective = 3
};

constexpr int noninjectiveMaxPreimages = 32;

class Evolvent
{
protected:
  int mDimension;
  int mTightness;

  std::vector<double> mRho;
  std::vector<double> mShiftScalars;
  std::vector<double> p2;

  void TransformToStandardCube(const double *y, double *z);
  void TransformToSearchDomain(const double *y, double *z);

  bool mIsInitialized;
private:
  MapType mMapType;
  int mMapKey;

public:
  Evolvent();
  Evolvent(int dimension, int tightness, const double* lb, const double* ub, MapType type = Simple);
  ~Evolvent();

  virtual void GetImage(double x, double y[]);
  virtual int GetAllPreimages(const double* p, double xp[]);
};

}

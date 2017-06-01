// Copyright (c) 2017 So Takamoto
// Released under the MIT license
// http://opensource.org/licenses/mit-license.php

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <type_traits>

#include <gmp.h>
#include <mpfr.h>

constexpr int mpfr_bit = 512;

using mpfr_t_arr2p = std::add_pointer<std::remove_extent<mpfr_t>::type>::type;
class mpfr_t_wrapper{
  friend mpfr_t_wrapper;
  mpfr_t_arr2p data;
  private:
  mpfr_t_wrapper(const mpfr_t_wrapper&);
  mpfr_t_wrapper& operator=(const mpfr_t_wrapper&);
  public:
  mpfr_t_wrapper(){
    data = new mpfr_t;
    mpfr_init2(data, mpfr_bit);
  }
  mpfr_t_wrapper(mpfr_t_wrapper&& r) noexcept{
    data = r.data;
    r.data = nullptr;
  }
  ~mpfr_t_wrapper(){
    if(data == nullptr) return;
    mpfr_clear(data);
    delete[] data;
  }
  operator mpfr_t_arr2p (){ return data;}
  operator const mpfr_t_arr2p () const{ return data;}
  mpfr_t_arr2p get(){ return data;}
  const mpfr_t_arr2p get() const{ return data;}
  mpfr_t_arr2p operator-> (){ return data;}
  const mpfr_t_arr2p operator-> () const{ return data;}
};

struct mass_struct{
  std::vector<mpfr_t_wrapper> m, l, k; //mass, length, spring constant
  mass_struct(size_t _n) : m(_n), l(_n), k(_n){}
  private:
  mass_struct (const mass_struct&);
  mass_struct& operator=(const mass_struct&);
};

struct point2d_t{
  mpfr_t x[2];
  point2d_t(){
    mpfr_init2(x[0], mpfr_bit);
    mpfr_init2(x[1], mpfr_bit);
  }
  ~point2d_t(){
    mpfr_clear(x[0]);
    mpfr_clear(x[1]);
  }
};

void accel(std::vector<point2d_t>& pos, const mass_struct& ms, std::vector<point2d_t>& ans){
  //initialize
  for(size_t i = 0; i < ans.size(); i++){
    mpfr_set_d(ans[i].x[0], 0.0, MPFR_RNDN);
    mpfr_set_d(ans[i].x[1], 0.0, MPFR_RNDN);
  }
  for(size_t i = 1; i < ans.size(); i++){
    mpfr_t_wrapper tmp, dr, drsq, f;
    point2d_t dx;
    //dx = pos[i]-pos[i-1]
    mpfr_sub(dx.x[0], pos[i].x[0], pos[i-1].x[0], MPFR_RNDN);
    mpfr_sub(dx.x[1], pos[i].x[1], pos[i-1].x[1], MPFR_RNDN);

    //potential energy (spring)
    //dr = sqrt(dx*dx)-ms.l[i]
    mpfr_sqr(drsq, dx.x[0], MPFR_RNDN);
    mpfr_sqr(tmp, dx.x[1], MPFR_RNDN);
    mpfr_add(drsq, drsq, tmp, MPFR_RNDN);
    mpfr_sqrt(dr, drsq, MPFR_RNDN);
    mpfr_sub(f, dr, ms.l[i], MPFR_RNDN);
    //f = ms.k[i]*dr
    mpfr_mul(f, f, ms.k[i], MPFR_RNDN);
    mpfr_div(f, f, dr, MPFR_RNDN);

    point2d_t fx;
    mpfr_mul(fx.x[0], f, dx.x[0], MPFR_RNDN);
    mpfr_mul(fx.x[1], f, dx.x[1], MPFR_RNDN);
    if(i != 1){
      mpfr_add(ans[i-1].x[0], ans[i-1].x[0], fx.x[0], MPFR_RNDN);
      mpfr_add(ans[i-1].x[1], ans[i-1].x[1], fx.x[1], MPFR_RNDN);
    }
    mpfr_sub(ans[i].x[0], ans[i].x[0], fx.x[0], MPFR_RNDN);
    mpfr_sub(ans[i].x[1], ans[i].x[1], fx.x[1], MPFR_RNDN);

    //potential energy (gravity)
    mpfr_mul_d(tmp, ms.m[i], -1.0, MPFR_RNDN);
    mpfr_add(ans[i].x[1], ans[i].x[1], tmp, MPFR_RNDN);
  }
  //acc = f/m
  for(size_t i = 0; i < ans.size(); i++){
    mpfr_div(ans[i].x[0], ans[i].x[0], ms.m[i], MPFR_RNDN);
    mpfr_div(ans[i].x[1], ans[i].x[1], ms.m[i], MPFR_RNDN);
  }
}

void energy(const std::vector<point2d_t>& pos, const std::vector<point2d_t>& dpos, const mass_struct& ms, mpfr_t_wrapper& ans){
  mpfr_t_wrapper tmp, total, dr;
  point2d_t dx;
  mpfr_set_d(total, 0.0, MPFR_RNDN);
  for(size_t i = 1; i < pos.size(); i++){
    //dx = pos[i]-pos[i-1]
    mpfr_sub(dx.x[0], pos[i].x[0], pos[i-1].x[0], MPFR_RNDN);
    mpfr_sub(dx.x[1], pos[i].x[1], pos[i-1].x[1], MPFR_RNDN);

    //potential energy (spring)
    //dr = sqrt(dx*dx)-ms.l[i]
    mpfr_sqr(dr, dx.x[0], MPFR_RNDN);
    mpfr_sqr(tmp, dx.x[1], MPFR_RNDN);
    mpfr_add(dr, dr, tmp, MPFR_RNDN);
    mpfr_sqrt(dr, dr, MPFR_RNDN);
    mpfr_sub(dr, dr, ms.l[i], MPFR_RNDN);
    //total += 0.5*ms.k[i]*dr*dr
    mpfr_sqr(dr, dr, MPFR_RNDN);
    mpfr_mul(dr, dr, ms.k[i], MPFR_RNDN);
    mpfr_mul_d(dr, dr, 0.5, MPFR_RNDN);
    mpfr_add(total, total, dr, MPFR_RNDN);

    //potential energy(gravity)
    //total += 1.0*pos[i].x[1]
    mpfr_mul(tmp, pos[i].x[1], ms.m[i], MPFR_RNDN);
    mpfr_mul_d(tmp, tmp, 1.0, MPFR_RNDN);
    mpfr_add(total, total, tmp, MPFR_RNDN);

    //kinetic energy
    //total += 0.5*ms.m[i]*dpos*dpos
    mpfr_sqr(dr, dpos[i].x[0], MPFR_RNDN);
    mpfr_sqr(tmp, dpos[i].x[1], MPFR_RNDN);
    mpfr_add(dr, dr, tmp, MPFR_RNDN);
    mpfr_mul(dr, dr, ms.m[i], MPFR_RNDN);
    mpfr_mul_d(dr, dr, 0.5, MPFR_RNDN);
    mpfr_add(total, total, dr, MPFR_RNDN);
  }
  mpfr_set(ans, total, MPFR_RNDN);
}

void euler(std::vector<point2d_t>& pos, std::vector<point2d_t>& dpos, const mass_struct& ms, const mpfr_t_wrapper& dt){
  mpfr_t_wrapper tmp;
  const int npos = pos.size();
  std::vector<point2d_t> ddpos(npos);
  accel(pos, ms, ddpos);
  for(int i = 0; i < npos; i++){
    for(int d = 0; d < 2; d++){
      mpfr_mul(tmp, dpos[i].x[d], dt, MPFR_RNDN);
      mpfr_add(pos[i].x[d], pos[i].x[d], tmp, MPFR_RNDN);

      mpfr_mul(tmp, ddpos[i].x[d], dt, MPFR_RNDN);
      mpfr_add(dpos[i].x[d], dpos[i].x[d], tmp, MPFR_RNDN);
    }
  }
}

void velocity_verlet_init(std::vector<point2d_t>& pos, std::vector<point2d_t>& dpos, const mass_struct& ms, std::vector<point2d_t>& current_ddpos){
  accel(pos, ms, current_ddpos);
}

void velocity_verlet(std::vector<point2d_t>& pos, std::vector<point2d_t>& dpos, const mass_struct& ms, const mpfr_t_wrapper& dt, std::vector<point2d_t>& current_ddpos){
  mpfr_t_wrapper tmp, dt2, dtsq2;
  const int npos = pos.size();

  mpfr_mul_d(dt2, dt, 0.5, MPFR_RNDN);
  mpfr_mul(dtsq2, dt2, dt, MPFR_RNDN);

  for(int i = 0; i < npos; i++){
    for(int d = 0; d < 2; d++){
      mpfr_mul(tmp, dpos[i].x[d], dt, MPFR_RNDN);
      mpfr_add(pos[i].x[d], pos[i].x[d], tmp, MPFR_RNDN);
      mpfr_mul(tmp, dtsq2, current_ddpos[i].x[d], MPFR_RNDN);
      mpfr_add(pos[i].x[d], pos[i].x[d], tmp, MPFR_RNDN);

      mpfr_mul(tmp, current_ddpos[i].x[d], dt2, MPFR_RNDN);
      mpfr_add(dpos[i].x[d], dpos[i].x[d], tmp, MPFR_RNDN);
    }
  }
  accel(pos, ms, current_ddpos);
  for(int i = 0; i < npos; i++){
    for(int d = 0; d < 2; d++){
      mpfr_mul(tmp, current_ddpos[i].x[d], dt2, MPFR_RNDN);
      mpfr_add(dpos[i].x[d], dpos[i].x[d], tmp, MPFR_RNDN);
    }
  }
}

void make_init_pos(std::vector<point2d_t>& pos, std::vector<point2d_t>& dpos, const mass_struct& ms, const mpfr_t_wrapper& diff){
  mpfr_set_d(pos[0].x[0], 0.0, MPFR_RNDN);
  mpfr_set_d(pos[0].x[1], 0.0, MPFR_RNDN);
  mpfr_set_d(dpos[0].x[0], 0.0, MPFR_RNDN);
  mpfr_set_d(dpos[0].x[1], 0.0, MPFR_RNDN);
  for(size_t i = 1; i < pos.size(); i++){
    mpfr_set_d(pos[i].x[0], 1.5, MPFR_RNDN);
    mpfr_set_d(pos[i].x[1], 1.5, MPFR_RNDN);
    if(i == pos.size()-1){
      mpfr_add(pos[i].x[0], pos[i].x[0], diff, MPFR_RNDN);
      mpfr_add(pos[i].x[1], pos[i].x[1], diff, MPFR_RNDN);
    }
    mpfr_cos(pos[i].x[0], pos[i].x[0], MPFR_RNDN);
    mpfr_sin(pos[i].x[1], pos[i].x[1], MPFR_RNDN);
    mpfr_mul(pos[i].x[0], pos[i].x[0], ms.l[i], MPFR_RNDN);
    mpfr_mul(pos[i].x[1], pos[i].x[1], ms.l[i], MPFR_RNDN);
    mpfr_add(pos[i].x[0], pos[i].x[0], pos[i-1].x[0], MPFR_RNDN);
    mpfr_add(pos[i].x[1], pos[i].x[1], pos[i-1].x[1], MPFR_RNDN);
    mpfr_set_d(dpos[i].x[0], 0.0, MPFR_RNDN);
    mpfr_set_d(dpos[i].x[1], 0.0, MPFR_RNDN);
  }
}

void save(const std::vector<std::vector<point2d_t>>& posArr, const std::vector<std::vector<point2d_t>>& dposArr, const mass_struct& ms, const int iloop, const int outstep, const bool file_save){
  const int nsystem = posArr.size();
  const int npoint = posArr[0].size();

  //Save as paraview .vtk file format
  if(file_save){
    std::string outstr = "out/pendulum."+std::to_string(iloop/outstep)+".vtk";
    std::fstream fw(outstr.c_str(), std::ios::out);
    fw << "# vtk DataFile Version 2.0\nTest1\nASCII\nDATASET POLYDATA\nPOINTS " << npoint*nsystem << " float\n";
    for(int j = 0; j < nsystem; j++){
      for(int k = 0; k < npoint; k++){
        double x = mpfr_get_d(posArr[j][k].x[0], MPFR_RNDN);
        double y = mpfr_get_d(posArr[j][k].x[1], MPFR_RNDN);
        fw << std::showpoint << x << " " << y << " 0.0\n";
      }
    }
    fw << "LINES " << (npoint-1)*nsystem << " " << 3*(npoint-1)*nsystem << "\n";
    for(int j = 0; j < nsystem; j++){
      int ninit = npoint*j;
      for(int k = 0; k < npoint-1; k++){
        fw << "2 " << ninit+k << " " << ninit+k+1 << "\n";
      }
    }
    fw << "POINT_DATA " << npoint*nsystem << "\n";
    fw << "SCALARS PID float 1\n";
    fw << "LOOKUP_TABLE default\n";
    for(int j = 0; j < nsystem; j++){
      for(int k = 0; k < npoint; k++){
        fw << j << "\n";
      }
    }
    fw.close();
  }

  double pos_diff_d = 0.0;
  if(posArr.size() >= 2){
    mpfr_t_wrapper pos_diff;
    mpfr_sub(pos_diff, posArr[0][npoint-1].x[0], posArr[1][npoint-1].x[0], MPFR_RNDN);
    pos_diff_d = mpfr_get_d(pos_diff, MPFR_RNDN);
  }
  mpfr_t_wrapper energy_0;
  energy(posArr[0], dposArr[0], ms, energy_0);
  double energy_d = mpfr_get_d(energy_0, MPFR_RNDN);
  std::cout << "Loop[" << iloop << "] energy: " << energy_d << ", diff: " << pos_diff_d << std::endl;
}

int main(int argc, char** argv){
  const int nsystem = 2;
  const int npoint = 4;

  double l = 1.5;
  mass_struct ms(npoint);
  for(int i = 0; i < npoint; i++){
    mpfr_set_d(ms.m[i], l, MPFR_RNDN);
    mpfr_set_d(ms.l[i], 1.0, MPFR_RNDN);
    mpfr_set_d(ms.k[i], 5000.0, MPFR_RNDN);
    l -= 0.4;
  }
  mpfr_set_d(ms.l[3], 0.5, MPFR_RNDN);

  std::vector<std::vector<point2d_t>> posArr(nsystem);
  std::vector<std::vector<point2d_t>> dposArr(nsystem);
  std::vector<std::vector<point2d_t>> ddposArr(nsystem);

  mpfr_t_wrapper small_diff;
  mpfr_set_d(small_diff, 1.0, MPFR_RNDN);
  //small_diff = 1.e-100
  for(int i = 0; i < 100; i++){
    mpfr_mul_d(small_diff, small_diff, 0.1, MPFR_RNDN);
  }

  mpfr_t_wrapper diff;
  mpfr_set_d(diff, 0.0, MPFR_RNDN);
  for(int i = 0; i < nsystem; i++){
    posArr[i] = std::vector<point2d_t>(npoint);
    dposArr[i] = std::vector<point2d_t>(npoint);
    make_init_pos(posArr[i], dposArr[i], ms, diff);

    ddposArr[i] = std::vector<point2d_t>(npoint);
    velocity_verlet_init(posArr[i], dposArr[i], ms, ddposArr[i]);

    mpfr_add(diff, diff, small_diff, MPFR_RNDN);
  }

  mpfr_t_wrapper timestep;
  mpfr_set_d(timestep, 0.0001, MPFR_RNDN);

  const bool file_save = false;

  const int outstep = 1000;
  const int nloop = 4000000;
  save(posArr, dposArr, ms, 0, outstep, file_save);
  for(int i = 0; i < nloop; i++){
    for(int j = 0; j < nsystem; j++){
      velocity_verlet(posArr[j], dposArr[j], ms, timestep, ddposArr[j]);
    }
    const int iloop = i+1;
    if(iloop%outstep == 0){
      save(posArr, dposArr, ms, iloop, outstep, file_save);
    }
  }

  /*
  //reverse dynamics calculation
  for(int j = 0; j < nsystem; j++){
    for(int k = 0; k < npoint; k++){
      mpfr_mul_d(dposArr[j][k].x[0], dposArr[j][k].x[0], -1.0, MPFR_RNDN);
      mpfr_mul_d(dposArr[j][k].x[1], dposArr[j][k].x[1], -1.0, MPFR_RNDN);
    }
    velocity_verlet_init(posArr[j], dposArr[j], ms, ddposArr[j]);
  }

  for(int i = nloop; i < 2*nloop; i++){
    for(int j = 0; j < nsystem; j++){
      velocity_verlet(posArr[j], dposArr[j], ms, timestep, ddposArr[j]);
    }
    const int iloop = i+1;
    if(iloop%outstep == 0){
      save(posArr, dposArr, ms, iloop, outstep, file_save);
    }
  }
  */

  return 0;
}


/*
 * soil_models.cc
 *
 *  Created on: Aug 11, 2016
 *      Author: jb
 */


#include <algorithm>
#include <cmath>

#include "fadbad.h"
#include "badiff.h"
#include "fadiff.h"

#include "system/asserts.hh"
#include "flow/soil_models.hh"


template <class Model>
void SoilModelImplBase<Model>::reset(SoilData data)
{
    model_.reset_(data);
}

template <class Model>
double SoilModelImplBase<Model>::conductivity( const double &p_head) const
{
    return model_.conductivity_(p_head);
}

template <class Model>
auto SoilModelImplBase<Model>::conductivity_diff(const DiffDouble &p_head)->DiffDouble const
{
    return model_.conductivity_(p_head);
}

template <class Model>
double SoilModelImplBase<Model>::water_content( const double &p_head) const
{
    return model_.water_content_(p_head);
}

template <class Model>
auto SoilModelImplBase<Model>::water_content_diff(const DiffDouble &p_head)->DiffDouble const
{
    return model_.water_content_(p_head);
}




namespace internal {

VanGenuchten::VanGenuchten()
: Bpar(0.5), Ppar(2), K_lower_limit(1.0E-20)
{}

void VanGenuchten::reset_(SoilData soil)
{
    // check soil parameters
    ASSERT_LT_DBG( soil.cut_fraction, 1.0);
    ASSERT_GT_DBG( soil.cut_fraction, 0.0);
    soil_param_ = soil;

    m = 1-1/soil_param_.n;
    // soil_param_.cut_fraction = (Qs -Qr)/(Qnc - Qr)
    Qs_nc = (soil_param_.Qs - soil_param_.Qr)/soil_param_.cut_fraction  + soil_param_.Qr;

    // conductivity internal scaling
    FFQr = 1.0;   // pow(1 - pow(Qeer,1/m),m);
    double Qs_relative = soil_param_.cut_fraction;
    FFQs = pow(1 - pow(Qs_relative, 1/m),m);


    Hs = Q_rel_inv(soil_param_.cut_fraction);
    //std::cout << "Hs : " << Hs << " qs: " << soil_param_.Qs << " qsnc: " << Qs_nc << std::endl;

}


template <class T>
T VanGenuchten::Q_rel(const T &h) const
{
    return  pow( 1 + pow(-soil_param_.alpha * h, soil_param_.n), -m);
}

template <class T>
T VanGenuchten::Q_rel_inv(const T &q) const
{
    return  -pow( pow( q, -1/m ) -1, 1/soil_param_.n)/soil_param_.alpha;
}



// pri generovani grafu jsem zjistil ze originalni vzorecek pro
// vodivost pres theta je numericky nestabilni pro tlaky blizke
// nule, zejmena pro velka n
// tento vzorec je stabilni:
//
/// k(h)=t(h)**(0.5)* (1- ((h)**n/(1+(h)**n)) **m)**2

template <class T>
T VanGenuchten::conductivity_(const T& h) const
{
    T Kr, Q_unscaled, Q_cut_unscaled, FFQ;

      if (h < Hs) {
            Q_unscaled = Q_rel(h);
            Q_cut_unscaled = Q_unscaled / soil_param_.cut_fraction;
            FFQ = pow(1 - pow(Q_unscaled,1/m),m);
            Kr = soil_param_.Ks * pow(Q_cut_unscaled,Bpar)*pow((FFQr - FFQ)/(FFQr - FFQs),Ppar);
    }
    else Kr = soil_param_.Ks;

    if (Kr < K_lower_limit) return K_lower_limit;
    else return Kr;
}

//template SoilModelBase::DiffDouble VanGenuchten::conductivity_<SoilModelBase::DiffDouble>(const SoilModelBase::DiffDouble &h) const;
//template double VanGenuchten::conductivity_<double>(const double &h) const;


template <class T>
T VanGenuchten::water_content_(const T& h) const
{
    if (h < Hs) return soil_param_.Qr + (Qs_nc - soil_param_.Qr) *Q_rel(h);
    else return soil_param_.Qs;
}

//template SoilModelBase::DiffDouble VanGenuchten::water_content_<SoilModelBase::DiffDouble>(const SoilModelBase::DiffDouble &h) const;
//template double VanGenuchten::water_content_<double>(const double &h) const;

Irmay::Irmay()
: VanGenuchten()
{
    Ppar=3;
}


template <class T>
T Irmay::conductivity_(const T& h) const
{
    T Kr;

    if (h < Hs) {
        T Q = this->Q_rel(h);
        Kr = soil_param_.Ks * pow(Q, Ppar);
    }
    else Kr = soil_param_.Ks;

    if (Kr < K_lower_limit) return K_lower_limit;
    else return Kr;
}

//template SoilModelBase::DiffDouble Irmay::conductivity_<SoilModelBase::DiffDouble>(const SoilModelBase::DiffDouble &h) const;
//template double Irmay::conductivity_<double>(const double &h) const;

} //close namespace internal



template class SoilModelImplBase<internal::VanGenuchten>;
template class SoilModelImplBase<internal::Irmay>;



/* 
 * File:   hydro_functions.hh
 * Author: jb
 *
 * Created on May 17, 2010, 1:52 PM
 */

#ifndef _HYDRO_FUNCTIONS_HH
#define	_HYDRO_FUNCTIONS_HH

#include <algorithm>
#include <cmath>

#include "base/parameter_handler.h"
using namespace dealii;


//class SoilModel {
//
//}

// todo: 
// objekt por SoilModel - jedna sada primarnich i pomocnych parametru
// konzistentni hydrologicke funkce
//
// potreboval bych mit vice modelu, z orezavanim i bez nej,
// orezavani zpusobuje velke residuum v miste nespojitosti derivace saturace

// problem - jak pouzit funktory, FADBAD, aproximace a 
// inverzi (inverze tlak v zav. na vlhkosti)
//
// podrobnejsi vyzkum ohledne hydrologickych funkci
// vodivost vyplyva s kapilarnich modelu (jak presne?)
// co vlhkost?
// pri generovani grafu jsem zjistil ze originalni vzorecek pro
// vodivost pres theta je numericky nestabilni pro tlaky blizke
// nule, zejmena pro velka n
// tento vzorec je stabilni:
//
/// k(h)=t(h)**(0.5)* (1- ((h)**n/(1+(h)**n)) **m)**2

class HydrologyModel {
public:
    HydrologyModel( ParameterHandler &prm );

    template <class T>
    T FK(const T &h) const;

    template <class T>
    T FQ(const T &h) const;

    /// Maximum point of capacity.
    inline double cap_arg_max() const { return arg_cap_max; }

    // function of normalized capacity (with maximum 1.0)
    inline double lambda(double cap) const {

        if (cap < 1.0/8) return 1;
        if (cap > 1.0/2) return 0.5;
        cap = ( cap - 1.0/8 ) /  (1.0/2 - 1.0/8);
        double lmb= //1-4*square(h/arg_cap_max);
                    0.5+0.5/
                            (1.0+square(
                                       1.0/(1.0- square(cap)) - 1.0
                                      )
                             );
        return lmb;

    }

private:
    inline double square(const double x) const
        {return x*x;}

    template <class T>
    inline T Q_non_cut(const T &h) const
    {
        if (h > 0.0) return Qs;
        else return Qr_nc + (Qs_nc - Qr_nc) * pow( 1 + pow(-alpha * h, n), -m);
    }

    template <class T>
    inline T Q_non_cut_inv(const T &q) const
    {
        return  -pow( pow( (q - Qr_nc) / (Qs_nc - Qr_nc), -1/m ) -1, 1/n)/alpha;
    }


    // input parameters

    double  n;              // power parameter
    double  alpha;           // pressure head scaling
    double  Qr;             // residual water content
    double  Qs;             // saturated water content
    double  cut_fraction;   // fraction of Qs where to cut and rescale
                        // water content and conductivity function
    double  Ks;             // saturated conductivity
    double  Kk;             // conductivity at the cut point

    // conductivity parameters
    const double Bpar;
    const double Ppar;
    const double K_lower_limit;

    // aux values
    double m;
    double arg_cap_max, cap_max_; // position of capacity maximum
    double Qs_nc;       // saturated value of continuation of cut water content function
    double Qr_nc;       // residual value of continuation of cut water content function

    double FFQr, FFQs;
    double Hs;

};

HydrologyModel:: HydrologyModel( ParameterHandler &prm )
: Bpar(0.5), Ppar(2), K_lower_limit(1.0E-20)
{
    n = prm.get_double("hydro_n");
    alpha = prm.get_double("hydro_alfa");

    cut_fraction = prm.get_double("hydro_cut");
    Qr = prm.get_double("hydro_qr");
    Qs = prm.get_double("hydro_qs");
    Ks = prm.get_double("hydro_ks");
    //Kk = prm.get_double("hydro_kk");


    m = 1-1/n;
    arg_cap_max= -alpha * pow(m, 1.0/n);

    Qs_nc = Qs / cut_fraction;
    Qr_nc = Qr; // no cut on residual part

    // conductivity internal scaling
    FFQr = 1.0;   // pow(1 - pow(Qeer,1/m),m);
    double Qs_unscaled = (Qs - Qr_nc) / (Qs_nc - Qr_nc );
    FFQs = pow(1 - pow(Qs_unscaled,1/m),m);


    Hs = Q_non_cut_inv(Qs);


}


template <class T>
T HydrologyModel::FK(const T& h) const
{
    T Kr,Q, Q_unscaled, Q_cut_unscaled, FFQ;

      if (h < Hs) {
            Q = FQ(h);
            Q_unscaled = (Q - Qr_nc) / (Qs_nc - Qr_nc);
            Q_cut_unscaled = (Q - Qr) / (Qs - Qr);

            FFQ = pow(1 - pow(Q_unscaled,1/m),m);

            Kr = Ks * pow(Q_cut_unscaled,Bpar)*pow((FFQr - FFQ)/(FFQr - FFQs),Ppar);
    }
    else Kr = Ks;

    if (Kr < K_lower_limit) return K_lower_limit;
    else return Kr;
}


template <class T>
T HydrologyModel::FQ(const T& h) const
{
    if (h < 0.0) return Q_non_cut(h);
    else return Qs;
}

/*

//FK-------------------------------------------------------------
template <class T>
class FK	//FK - hydraulic conductivity function
{
private:
// function parameters
static const double	Bpar;   // = 0.5;
static const double	PPar;   // = 2;

HydrologyParams par;

// auxiliary parameters
double Qa,Qm,Qk,m,Hr,Hk,Hs,C1Qee,C2Qee,Qeer,Qees,Qeek;
public:
    FK(const HydrologyParams &par);
    T operator()(const T&  h);
};

template <class T> const double FK<T>::Bpar =0.5;
template <class T> const double FK<T>::PPar =2;

template <class T>
FK<T> :: FK(const HydrologyParams &par)
: par(par)
{

	m = 1-1/par.n;

        Qa=par.Qr;
        Qm=par.Qs / par.cut_fraction;
        //Qk=0.99 * par.Qs;
        Qk=par.Qs;
        C1Qee = 1/(Qm - Qa);
        C2Qee = -Qa/(Qm - Qa);
//    if (par.cut_fraction >0.9999) {
//        Hr=-1000;
//        Hk=Hs=0;
//        return;
//    }
                                                        // fraction =1
        Qees = C1Qee * par.Qs + C2Qee;                  // 1
        Qeek = min(C1Qee * Qk + C2Qee, Qees);           // 1
        Qeer = min(C1Qee * par.Qr + C2Qee, Qeek);       // 0

        Hr = -1/par.alfa*pow(pow(Qeer,-1/m)-1,1/par.n); //
        Hs = -1/par.alfa*pow(pow(Qees,-1/m)-1,1/par.n);
        Hk = -1/par.alfa*pow(pow(Qeek,-1/m)-1,1/par.n);

        cout << "K: Hr, Hs, Hk:" << Hr << " " << Hs << " " << Hk << endl;
}

template <class T>
T FK<T> :: operator()(const T& h)
{
    T Kr,FFQr,FFQ,FFQk,Qee,Qe,Qek,C1Qe,C2Qe,Q;

    if (h < Hr) return par.Ks*(1E-200);      // numericaly zero Kr
    else 
      if (h < Hk) {
            Q = Qa + (Qm - Qa)*pow((1 + pow(-par.alfa*h,par.n)),-m);
            Qee = C1Qee*Q + C2Qee;
            FFQr = pow(1 - pow(Qeer,1/m),m);
            FFQ = pow(1 - pow(Qee,1/m),m);
            FFQk = pow(1 - pow(Qeek,1/m),m);
            C1Qe = 1/(par.Qs - par.Qr);
            C2Qe = -par.Qr/(par.Qs - par.Qr);
            Qe = C1Qe*Q + C2Qe;
            Qek = C1Qe*Qk + C2Qe;
            Kr = pow(Qe/Qek,Bpar)*pow((FFQr - FFQ)/(FFQr - FFQk),PPar) * par.Kk/par.Ks;
            return ( max<T>(par.Ks*Kr,par.Ks*(1E-10)) );
    }
    else if(h <= Hs)
    {
         Kr = (1-par.Kk/par.Ks)/(Hs-Hk)*(h-Hs) + 1;
         return ( par.Ks*Kr );
    }
    else return par.Ks;
}


//FQ--------------------------------------------------------------
template <class T>
class FQ	//FQ - water saturation function
{
HydrologyParams par;
// auxiliary parameters
double  m, Qeer,Qees,Hr,Hs, Qa, Qm;
public:
    FQ(const HydrologyParams &par);
    T operator()(const T&  h);
};

template <class T>
FQ<T>::FQ(const HydrologyParams &par)
: par(par)
{
    m = 1 - 1 / par.n;

//    if (par.cut_fraction >0.9999) {
//        Hr=-1000;
//        Hs=0;
//        return;
//    }
    Qm=par.Qs / par.cut_fraction;
    Qa=par.Qr;

    Qeer = (par.Qr - Qa)/ (Qm - Qa);
    Qees = (par.Qs - Qa)/ (Qm - Qa);
    Hr = -1 / par.alfa * pow( pow(Qeer,-1/m)-1, 1/par.n);
    Hs = -1 / par.alfa * pow( pow(Qees,-1/m)-1, 1/par.n);

    cout << "Q: Hr, Hs:" << Hr << " " << Hs << " " << endl;
}


template <class T>
T FQ<T>::operator()(const T& h)
{
    T  Qee;

    if (h < Hr) return par.Qr;
    else if (h < Hs) {
        Qee = pow( 1 + pow(-par.alfa * h, par.n), -m);
        return Qa + (Qm - Qa) * Qee;
    }
    else return par.Qs;

}

//FC--------------------------------------------------------------
template <class T>
class FC	//FC - water capacity function FQ derivative
{
HydrologyParams par;
// auxiliary parameters
double  m, C1Qee,C2Qee,Qeer,Qees,Hr,Hs, Qa, Qm;
public:
    FC(const HydrologyParams &par);
    T operator()(const T&  h);
};

template <class T>
FC<T>::FC(const HydrologyParams &par)
: par(par)
{
    m = 1 - 1 / par.n;
    Qm=par.Qs / par.cut_fraction;
    Qa=par.Qr;

    C1Qee = 1/(Qm - Qa);
    C2Qee = -Qa/(Qm - Qa);
    Qeer = C1Qee * par.Qr + C2Qee;
    Qees = C1Qee * par.Qs + C2Qee;
    Hr = -1 / par.alfa * pow( pow(Qeer,-1/m)-1, 1/par.n);
    Hs = -1 / par.alfa * pow( pow(Qees,-1/m)-1, 1/par.n);
}


template <class T>
T FC<T>::operator()(const T& h)
{
    T  C1;

    if (h <= Hr)
        return 0.0;
    else if ( h < Hs ) {
        C1=pow( 1 + pow( -par.alfa * h, par.n), -m-1);
        return (Qm - Qa) * m * par.n * pow(par.alfa, par.n)*pow(-h, par.n-1)*C1;
    }
    else
        return 0.0;
}
/*

          ! ************************************************************************
! FH - inverse water saturation function
! **************************************
real pure function FH_4(h,Matr)
implicit none
integer,INTENT(IN) :: Matr
real, INTENT(IN) :: h
  FH_4=FH_8(DBLE(h),Matr)
end function FH_4

double precision pure function FH_8(Qe,Matr)
implicit none
integer, INTENT(IN) :: Matr
double precision,INTENT(IN) :: Qe
double precision :: n,m,Qr,Qs,Qa,Qm,Alfa,Q,Qee

  Qr=par(1,Matr)
  Qs=par(2,Matr)
  Qa=par(3,Matr)
  Qm=par(4,Matr)
  Alfa=par(5,Matr)
  n=par(6,Matr)
  m=1-1/n
  Q=Qr+(Qs-Qr)*Qe
  Qee=dmax1((Q-Qa)/(Qm-Qa),1d-3)
  Qee=dmin1(Qee,1-1d-15)
  FH_8=-1/Alfa*(Qee**(-1/m)-1)**(1/n)
end function FH_8
*/

// Conductivity for analytical solution

class HydroModel_analytical
{
public:
    HydroModel_analytical(ParameterHandler &prm) {};

    /// Maximum point of capacity. (computed numericaly from eq. atan(x)*2x=1 )
    inline double cap_arg_max() const { return -0.765378926665788882857; }

    template <class T>
    T FK(const T &h) const
    {
        if (h>=0.0) return ( 2.0 );
        return ( 2.0 / (1+h*h) );
    }

    template <class T>
    T FQ(const T &h) const
    {
        static T pi_half =std::atan(1.0)*2;
        if (h>=0.0) return ( 2.0*pi_half*pi_half );
        T a_tan = atan(h);
        return ( 2*pi_half*pi_half - a_tan*a_tan );
    }
private:
    double arg_cap_max;
};


class HydroModel_linear
{
    HydroModel_linear(ParameterHandler &prm) {};

public:
    template <class T>
    T FK(const T &h) const {
        return 1.0;
    }

    template <class T>
    T FQ(const T &h) const {
        return h;
    }
};


#endif	/* _HYDRO_FUNCTIONS_HH */


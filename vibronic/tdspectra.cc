#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include "constants.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

using namespace std;
using namespace psi;

// From Eigen documentation:
//  class Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >
//  Some notes:
//   Dynamic (predefined typedef): Dynamic size
//   RowMajor (option): Storage order is row majo (row by row)
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> EigenVector;
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixC;
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> EigenVectorC;


void tdspectra(double DE, Vector &G1_,Vector &G2_,Matrix &GM_,Vector &GK_,
               double* dip0, double* dipm0, double** dipq, double** dipmq,
               Options& options)
{

    // Get options
    double temp               = options.get_double("TEMPERATURE");
    std::string opdip         = options.get_str("DIPOLE_APPROX");
    std::string optrans       = options.get_str("TRANSITION");
    std::string optrans1      = options.get_str("CHIRAL");
    std::string opmod         = options.get_str("MODEL_PES");
    int corr_npoints          = options.get_int("CORR_NPOINTS");
    double corr_time          = options.get_double("CORR_TIME");
    double inh_brd_hwhm       = options.get_double("BROAD_HWHM");
    std::string inh_brd_type  = options.get_str("BROAD_TYPE");
    // Convert to AU 
    inh_brd_hwhm /= autoev;
    corr_time    /= autofs;


    // Transform psi::Matrix/Vector into Eigen
    int Nvib = GM_.ncol();
    EigenMatrix GM(Nvib,Nvib);
    EigenVector GK(Nvib), G1(Nvib), G2(Nvib);
    for (int i=0; i<Nvib; i++) {
        GK(i) = GK_(i);
        G1(i) = G1_(i);
        G2(i) = G2_(i);
        for (int j=0; j<Nvib; j++) {
            GM(i,j) = GM_(i,j);
        }
    }

    // Intermediate files for debuging
    ofstream cc51("cc.51");
    ofstream cc53("cc.53");
    ofstream cc54("cc.54");
    ofstream cc55("cc.55");
    ofstream cc56("cc.56");
    ofstream cc57("cc.57");
    
    // Initialize variables
    // Corr function (using fftw classes)
    double* t         = new double[2*corr_npoints];
    fftw_complex corrtime[2*corr_npoints];
    //
    double dtime=2.*corr_time/(double)(2*corr_npoints-1);
    double boltz = temp*kboltz/HARTtoJ;
    complex<double> Im(0, 1);
    complex<double> Re(1, 0);
    // Inh. broadening
    double aexp;
    int nbroad;
    if (inh_brd_type == "GAUSSIAN") {
        aexp = inh_brd_hwhm/sqrt(2.*log(2.));
        aexp = pow(aexp,2)/2.;
        nbroad = 2;
    }
    else if (inh_brd_type == "LORENTZIAN") {
        aexp = inh_brd_hwhm;
        nbroad = 2;
    }
    
    // Cicle around time
    complex<double> det_prev = Re;
    for (int it=0; it<corr_npoints; it++) {
        double time = (double)(it)*dtime+dtime/2.;
        complex<double> tauf = time;
        complex<double> taui;
        if (temp != 0.) {
            taui =-time-Im/boltz;
        }
        
        EigenMatrixC C_TD = EigenMatrixC::Zero(Nvib,Nvib);
        EigenMatrixC D_TD = EigenMatrixC::Zero(Nvib,Nvib);
        EigenVectorC asmfTD(Nvib);
        EigenVectorC asmiTD(Nvib);
        EigenVectorC csmfTD(Nvib);
        EigenVectorC csmiTD(Nvib);
        EigenVectorC dsmfTD(Nvib);
        EigenVectorC dsmiTD(Nvib);
        for (int i=0; i<Nvib; i++)  {
            if (temp == 0.) {
                    asmfTD(i) = G2(i)/sin(G2(i)*tauf);
                    asmiTD(i) = -2*G1(i)/(sin(G1(i)*time) + Im*cos(G1(i)*time));
                    csmfTD(i) = G2(i)                                                
                                 //*coth(Im*G2(i)*tauf/2.)                           
                                 *(exp(Im*G2(i)*tauf/2.)+exp(-Im*G2(i)*tauf/2.)) 
                                 /(exp(Im*G2(i)*tauf/2.)-exp(-Im*G2(i)*tauf/2.));
                    csmiTD(i) = G1(i);
                    dsmfTD(i) = G2(i)                            
                                 //*tanh(Im*G2(i)*tauf/2.) 
                                 *(exp(Im*G2(i)*tauf/2.)-exp(-Im*G2(i)*tauf/2.)) 
                                 /(exp(Im*G2(i)*tauf/2.)+exp(-Im*G2(i)*tauf/2.));
                    dsmiTD(i) = G1(i);

                    C_TD(i,i) = csmfTD(i);
                    D_TD(i,i) = dsmfTD(i);
            }
            else { // finite temperature
                asmfTD(i) = G2(i)/sin(G2(i)*tauf);
                // asmiTD and Z(part function) together
                asmiTD(i) = -2.*G1(i)/                                         
                             ( sin(G1(i)*time)/(1.-1./cos(Im*G1(i)/boltz)) + 
                               cos(G1(i)*time)/( 1./tan(Im*G1(i)/boltz)    - 
                               1./sin(Im*G1(i)/boltz)));
                csmfTD(i) = G2(i)                                           
                             //*coth(Im*G2(i)*tauf/2.)                        
                             *(exp(Im*G2(i)*tauf/2.)+exp(-Im*G2(i)*tauf/2.)) 
                             /(exp(Im*G2(i)*tauf/2.)-exp(-Im*G2(i)*tauf/2.));
//                              //Alternative formulation of coth
//                              *(1.+exp(-Im*G2(i)*tauf))                       
//                              /(1.-exp(-Im*G2(i)*tauf));
//                     cout << G2(i) << "*exp(Im*" << G2(i) << tauf/2. << ")+"
//                                   << "exp(-Im*" << G2(i) << tauf/2. << "))"
//                                   << "*exp(Im*" << G2(i) << tauf/2. << ")-"
//                                   << "exp(-Im*" << G2(i) << tauf/2. << "))";
//                     cout <<endl<<endl;
                csmiTD(i) = G1(i)      
                             //*coth(Im*G1(i)*taui/2.)
                             *(exp(Im*G1(i)*taui/2.)+exp(-Im*G1(i)*taui/2.)) 
                             /(exp(Im*G1(i)*taui/2.)-exp(-Im*G1(i)*taui/2.));
                             //Alternative formulation of coth
//                              *(1.+exp(-Im*G1(i)*taui))                       
//                              /(1.-exp(-Im*G1(i)*taui));
                dsmfTD(i) = G2(i)                                            
                             //*tanh(Im*G2(i)*tauf/2.)
                             *(exp(Im*G2(i)*tauf/2.)-exp(-Im*G2(i)*tauf/2.)) 
                             /(exp(Im*G2(i)*tauf/2.)+exp(-Im*G2(i)*tauf/2.));
                dsmiTD(i) = G1(i) 
                             //*tanh(Im*G1(i)*taui/2.)
                             *(exp(Im*G1(i)*taui/2.)-exp(-Im*G1(i)*taui/2.)) 
                             /(exp(Im*G1(i)*taui/2.)+exp(-Im*G1(i)*taui/2.));
                             //Alternative formulation of tanh
//                              *(1.-exp(-Im*G1(i)*taui))                       
//                              /(1.+exp(-Im*G1(i)*taui));
                             
                C_TD(i,i) = csmfTD(i);
                D_TD(i,i) = dsmfTD(i);
            }
        }
        // Complete the calculation of C/D matrices taking advantage of symmetry
        for (int i=0; i<Nvib; i++) {
            for (int j=0; j<=i; j++) {
                complex<double> auxC = (0.,0.);
                complex<double> auxD = (0.,0.);
                for (int k=0; k<Nvib; k++) {
                    auxC += csmiTD(k)*GM(k,j)*GM(k,i);
                    auxD += csmiTD(k)*GM(k,j)*GM(k,i);
                }
                C_TD(i,j) += auxC;
                D_TD(i,j) += auxD;
                // Are symmetric
                C_TD(j,i) = C_TD(i,j);
                D_TD(j,i) = D_TD(i,j);
                

            }
        }
        
        // Determinant
        complex<double> det = 1./C_TD.determinant()/D_TD.determinant();
        for (int i=0; i<Nvib; i++) {
            det *= asmiTD(i)*asmfTD(i);
            det *= -1.;
        }
        det = sqrt(det);
        
        double a_test = abs(det.real()-det_prev.real()) +
                        abs(det.imag()-det_prev.imag());
        double b_test = abs(det.real()+det_prev.real()) +
                        abs(det.imag()+det_prev.imag());
        if (b_test < a_test) {
            det = -det;
        }
        det_prev = det;
        
        cc51 << time*autofs<<"  " << det.real()<<"  " << det.imag()<<endl;
        
        // COMPUTE EXPONENTIAL
        complex<double> corr = (0.,0.);

        // K^tdfK
        for (int i=0; i<Nvib; i++) {
            corr -= GK(i)*dsmiTD(i)*GK(i);
        }
        // K^tdfJ
        EigenVectorC Vaux2 = EigenVectorC::Zero(Nvib);
        for (int i=0; i<Nvib; i++) {
            for (int j=0; j<Nvib; j++) {
                Vaux2(i) += GK(j)*dsmiTD(j)*GM(j,i);
            }
        }
        // J^tdiK
        EigenVectorC Vaux = EigenVectorC::Zero(Nvib);
        for (int i=0; i<Nvib; i++) {
            for (int j=0; j<Nvib; j++) {
                Vaux(i) += GM(j,i)*dsmiTD(j)*GK(j);
            }
        }
        // D^-1
        D_TD = D_TD.inverse();
        
        // Gather all terms into corr
        for (int i=0; i<Nvib; i++) {
            for (int j=0; j<Nvib; j++) {
                corr += Vaux2(i)*D_TD(i,j)*Vaux(j);
            }
        }
                
        // Track exponential and its argument
        cc53 << time*autofs<<"  " << exp(corr).real()<<"  " << exp(corr).imag()<<endl;

        //Also Ead and broad functions to reconstruct the whole correlation function
        cc54 << time*autofs<<"  " << exp(-Im*DE*time).real()<<"  " << exp(-Im*DE*time).imag()<<endl;
        cc55 <<  time*autofs<<"  " << exp(-aexp*pow(time,nbroad))<<endl;
        
        //----------------------
        // COMPUTE DIPOLE TERM
        //----------------------
        // FC
        complex<double> dipTerm = dip0[0]*dipm0[0]+dip0[1]*dipm0[1]+dip0[2]*dipm0[2];
        if (opdip != "FC") {
            EigenVectorC D_HT = EigenVectorC::Zero(Nvib);
            //Recall that we already have:
            // D^-1   stored in D_TD
            // J^tdfK stored in Vaux
            for (int i=0; i<Nvib; i++) {
                for (int j=0; j<Nvib; j++) {
                    D_HT(i) += D_TD(i,j)*Vaux(j);
                }
            }
            // C^-1 (store in C_TD)
            C_TD = C_TD.inverse();
            // A_HT
            EigenMatrixC A_HT(Nvib,Nvib);
            for (int i=0; i<Nvib; i++) {
                for (int j=0; j<=i; j++) {
                    A_HT(i,j) = (D_TD(i,j)-C_TD(i,j))/2.+D_HT(i)*D_HT(j);
                    A_HT(j,i) = A_HT(i,j);
                }
            }
            // HT contributions
            for (int i=0; i<Nvib; i++) {
                // FC/HT
                dipTerm += (dip0[0]*dipmq[i][0] + dip0[1]*dipmq[i][1]+dip0[2]*dipmq[i][2])*D_HT(i) +
                           (dipm0[0]*dipq[i][0] + dipm0[1]*dipq[i][1]+dipm0[2]*dipq[i][2])*D_HT(i);
                // HT
                for (int j=0; j<Nvib; j++) {
                    dipTerm += (dipq[i][0]*dipmq[j][0]+dipq[i][1]*dipmq[j][1] +
                               dipq[i][2]*dipmq[j][2])*A_HT(i,j);
                }
            }
        }

        // Track dipole term
        cc56 << time*autofs<<"  " << dipTerm.real()<<"  " << dipTerm.imag()<<endl;
        
        
        //-----------------------------
        // COMPUTE SWITCH FUNCTION (GIBBS)
        //-----------------------------
        double damp = (pow(cos(PI/2.*time/corr_time),2));
        cc55 <<  time*autofs<<"  " << damp<<endl;
            

        //=======================
        // CORRELATION FUNCTION:
        //=======================
        corr = exp(corr)*exp(-Im*DE*time)*exp(-aexp*pow(time,nbroad))*det*dipTerm*damp;
        
        //-----------------------------
        // Save in memory
        //-----------------------------
        //   * Time is stored in fs
        //   * Xi(t) in atomic units
        // Store from -t to t uisng Hemiticity: Xi(t) = [Xi(-t)]*
        // Positive X-axis
        t[it+corr_npoints]        = (it*dtime+dtime/2.)*autofs;   
        corrtime[it+corr_npoints][0] = corr.real();     
        corrtime[it+corr_npoints][1] = corr.imag();
        // Negative X-axis (hermitian property)
        t[corr_npoints-it-1]         = (-it*dtime-dtime/2.)*autofs;      
        corrtime[corr_npoints-it-1][0]  = corr.real(); 
        corrtime[corr_npoints-it-1][1]  =-corr.imag();

    }
    cc51.close();
    cc53.close();
    cc54.close();
    cc55.close();
    cc56.close();
    cc57.close();
    
    ofstream tcorr("timecorr.dat");
    for (int i=0; i<2*corr_npoints; i++) {
        tcorr << t[i] <<" "<<corrtime[i][0]<<" "<<corrtime[i][1]<<endl;
    }
    tcorr.close();
    
    // FFTW
    fftw_plan plan;
    int N=2*corr_npoints;
    fftw_complex spec_fftw[N];
    plan = fftw_plan_dft_1d(N,corrtime,spec_fftw,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    // Convert from DFT to FFT
    complex<double>* spec = new complex<double>[N];
    for (int i=0; i<N; i++) {
        spec[i] = (spec_fftw[i][0],spec_fftw[i][1]);
    }
    double df = 1./(t[N-1]-t[0]);
    double t0 = t[0];
    for (int k; k<N/2; k++) {
        spec[k] *= dtime * exp(-Im * 2.*PI*df * t0 * (double)(k-1)) /autofs;
        spec[k] /= 2.*PI;
    }
    
    ofstream spc("spectrum.dat");
    for (int i=0; i<corr_npoints; i++) {
        spc << (double)i*df*fstoev <<" "<<spec[i].real()<<" "<<spec[i].real()<<endl;
    }
    spc.close();
    
    return;
}


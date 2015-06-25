#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstring>

using namespace psi;

Matrix rot_to_fit_rmsd(Matrix &geo1, Matrix &geo2)
{  
    int N = geo1.nrow();
    
    Matrix ar1(4,4);
    Matrix al2(4,4);
    Matrix aaq(4,4);   
    
    // Build Quaternion matrix (aaq)
    aaq.zero();
    for (int i=0; i<N; i++) {
        ar1(0,0)=0.;           
        ar1(1,0)= geo1(i,0);    
        ar1(2,0)= geo1(i,1);    
        ar1(3,0)= geo1(i,2);    
        ar1(0,1)=-geo1(i,0);   
        ar1(1,1)=0.;           
        ar1(2,1)=-geo1(i,2);   
        ar1(3,1)= geo1(i,1);    
        ar1(0,2)=-geo1(i,1);   
        ar1(1,2)= geo1(i,2);    
        ar1(2,2)=0.;           
        ar1(3,2)=-geo1(i,0);   
        ar1(0,3)=-geo1(i,2);   
        ar1(1,3)=-geo1(i,1);   
        ar1(2,3)= geo1(i,0);    
        ar1(3,3)=0.;           
 
        al2(0,0)=0.;           
        al2(1,0)= geo2(i,0);    
        al2(2,0)= geo2(i,1);    
        al2(3,0)= geo2(i,2);    
        al2(0,1)=-geo2(i,0);   
        al2(1,1)=0.;           
        al2(2,1)= geo2(i,2);    
        al2(3,1)=-geo2(i,1);   
        al2(0,2)=-geo2(i,1);   
        al2(1,2)=-geo2(i,2);   
        al2(2,2)=0.;           
        al2(3,2)= geo2(i,0);    
        al2(0,3)=-geo2(i,2);   
        al2(1,3)= geo2(i,1);    
        al2(2,3)=-geo2(i,0);   
        al2(3,3)=0.;           
       
        aaq.accumulate_product(&al2,&ar1);
    }
    aaq.scale(-1.);
    
    // Diagonalize Quaternion matrix
    Matrix eigenvec(4,4);
    Vector eigenval(4);
    aaq.diagonalize(&eigenvec,&eigenval);
    
    // Compute the rotation matrix
    ar1.zero();
    ar1.accumulate_product(&aaq,&eigenvec);
    aaq.zero();
    aaq.accumulate_product(&eigenvec,&ar1);
    
    Matrix rot(3,3);
    rot(0,0)=pow(eigenvec(0,3),2)+pow(eigenvec(1,3),2)-pow(eigenvec(2,3),2)-pow(eigenvec(3,3),2);
    rot(1,0)=2.*(eigenvec(1,3)*eigenvec(2,3)+eigenvec(0,3)*eigenvec(3,3));
    rot(2,0)=2.*(eigenvec(1,3)*eigenvec(3,3)-eigenvec(0,3)*eigenvec(2,3));
    rot(0,1)=2.*(eigenvec(1,3)*eigenvec(2,3)-eigenvec(0,3)*eigenvec(3,3));
    rot(1,1)=pow(eigenvec(0,3),2)-pow(eigenvec(1,3),2)+pow(eigenvec(2,3),2)-pow(eigenvec(3,3),2);
    rot(2,1)=2.*(eigenvec(2,3)*eigenvec(3,3)+eigenvec(0,3)*eigenvec(1,3));
    rot(0,2)=2.*(eigenvec(1,3)*eigenvec(3,3)+eigenvec(0,3)*eigenvec(2,3));
    rot(1,2)=2.*(eigenvec(2,3)*eigenvec(3,3)-eigenvec(0,3)*eigenvec(1,3));
    rot(2,2)=pow(eigenvec(0,3),2)-pow(eigenvec(1,3),2)-pow(eigenvec(2,3),2)+pow(eigenvec(3,3),2);

    // We also rotate state1 in this subroutine
    geo1.transform(rot);

    return rot;
        
}

void rotate3d_3Nmat(Matrix& T, Matrix& rot) 
{
    int Nvib = T.ncol();
    int Nat  = T.nrow()/3;
    Matrix A(3*Nat,Nvib);
    
    for (int i=0; i<3*Nat; i++) {
        for (int j=0; j<Nvib; j++) {
            // Row of R: i%3 
            int ii = i%3;
            // Col of T: j, but only partial: [(i-i%3 : i-i%3+2),j] 
            int k  = i - i%3;
            A(i,j) = rot(ii,0)*T(k+0,j) + 
                     rot(ii,1)*T(k+1,j) + 
                     rot(ii,2)*T(k+2,j);
        }
    }
    
    T.copy(A);
    
    return;
}


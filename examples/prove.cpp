#include "stdinclude.hpp"
#include "POD_basis.hpp"
#include "rw_restart.hpp"

int main() {

    complex<double> a(1.0, -2.0);
    complex<double> b(2.0, 4.0);
    complex<double> c(0.0, -7.0);
    complex<double> d(5.0, 3.0);

    MatrixXcd C(3,2);
    C(0,0) = a;
    C(0,1) = b;
    C(1,0) = c;
    C(1,1) = d;
    C(2,0) = 2.0*c;
    C(2,1) = 0.5*d;

    cout << "C: \n" << C << endl;
    cout << "Transpose of C: \n" << C.conjugate() << endl; 
    cout << "C*C: \n" << C.transpose()*C << endl;
    // cout <<  "z is : " << z << endl;
    // cout << "The complex logarithm of z is : " << log(z) << endl;

    // VectorXcd zv(1);
    // zv(0) = z;
    // cout << "The real part of z is : " << zv.real() << endl;

    // cout << "Ciao ciao" << endl;


     return 0;

    //int Nr = 4;
    //int Ns = 3;
    //MatrixXd snap = MatrixXd::Random(Nr, Ns);
    //cout << "snap random: \n" << snap << endl << endl;
//
    //
//
    //VectorXd mean = snap.rowwise().mean();
    //int kk;
    //int nt = 1.5;
    //double t = 0.0;
//
    //VectorXd K_pc(Ns), lam(Ns);
    //MatrixXd Coeffs(Ns, Ns);
//
    //// MatrixXd phi = POD_basis( Nr, snap, K_pc, lam, Coeffs);
    //// 
    //// cout << " Optimal basis: \n " << phi << endl << endl;
//
    //double t_init = 0.0;
    //double Dt = 0.5;
    //double t_star = 0.0;
//
    //// VectorXd coefs_interp = RBF_Coefs( Coeffs.transpose(), phi, Dt, t_star, t_init);
//
    //// VectorXd Rec = phi*coefs_interp + mean;
//
    //// cout << " Reconstruction: \n " << Rec << endl;
//
    //VectorXd Init_state = snap.col(0);
    //VectorXcd lam_dmd;
//
    //// VectorXd lam(Ns), K_pc(Ns);
    //// MatrixXd Coeffs(Ns, Ns);
    //// MatrixXd phi_POD = POD_basis( Nr, snap.leftCols(Ns-1), K_pc, lam, Coeffs, "DMD");
    //// cout << " Lambda POD : \t" << lam.transpose() << endl;
    //// cout << "Number of snaps : " << Ns << endl;
    //// cout << " Number of non zero modes : " << phi_POD.cols() << endl;
//
//
//
    //cout << "-------------------DMD----------------------" << endl;
//
//
//
//
    //MatrixXcd phi_dmd = DMD_basis( Nr, snap, lam_dmd);
//
    //cout << " DMD modes : \n" << phi_dmd << endl;
//
    //VectorXcd b = coefs_dmd( phi_dmd, Init_state );
//
    //cout << "Lambda DMD : \t" << lam_dmd.transpose() << endl;
    //cout << "coefs DMD : \t" << b.transpose() << endl;
    //cout << " Initial state : \n" << Init_state << endl << endl;
//
//
//
    //VectorXcd Rec_dmd = Rec_field_dmd ( Nr, phi_dmd, lam_dmd, b, t );
    //VectorXd Recon = Rec_dmd.real();
//
    //cout << "dmd result: \n " << Rec_dmd << endl << endl;
    //cout << " dmd result without function rec : \n" << phi_dmd*b << endl << endl;
//
    //// cout << " mean : \n" << mean << endl;
//
    //return 0;




}
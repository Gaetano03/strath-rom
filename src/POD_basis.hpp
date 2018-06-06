#ifndef POD_HPP
#define POD_HPP

#include "stdinclude.hpp"
#include "rw_restart.hpp"


using namespace smartuq;
using namespace surrogate;




MatrixXd POD_basis( int Nr, MatrixXd snap, VectorXd &K_pc, VectorXd &lam, MatrixXd &eig_vec, string flag = "POD"){


    int N_snap = snap.cols();

    // MatrixXd field(Nr, 1);
    MatrixXd R(N_snap, N_snap);
    MatrixXd phi_c(Nr,N_snap);
    // MatrixXd eig_vec(N_snap, N_snap);

    VectorXd mean(Nr);
 
    mean = snap.rowwise().mean();

    if ( flag == "POD" ){
        for (int i = 0; i < N_snap; i++)
            snap.col(i) -= mean;
    }

    R = snap.transpose()*snap;

    EigenSolver<MatrixXd> es(R); 
    lam = es.eigenvalues().real();
    eig_vec = es.eigenvectors().real();
    eig_sort(lam, eig_vec);


    // for ( int i = 0; i < N_snap ; i++){
    //     cout << lam(i) << "\t" << eig_vec.row(i) << endl; 
    // }

    // cout << "Eigenvalues:\n " << lam << endl << endl;
    // cout << "EIG_T*EIG : " << endl;
    // MatrixXd eigcheckort = eig_vec.transpose()*eig_vec;
    // for ( int i = 0; i < N_snap; i++ ){
    //     for ( int j = 0; j < N_snap; j++ ){
    //         if ( eigcheckort(i,j) < 1e-12 )
    //             eigcheckort(i,j) = 0.0;
    //     }

    // } 
    // cout << eigcheckort << endl;


    double sum = 0;

    for (int i = 0; i < N_snap; i++){
        sum += lam(i)/lam.sum();
        K_pc(i) = sum;
    }

    // cout << "Energy content : \n" << K_pc << endl << endl;
    
    double tol = lam(0)*1e-12;
    phi_c = snap*eig_vec;

    int count = 0;
    while ( count < lam.size() && abs(lam(count)) > tol)
            count++;

    MatrixXd phi(Nr,count);
    for ( int i = 0 ; i < count ; i++ ){
            phi.col(i) = phi_c.col(i)/sqrt(lam(i));
    
    }

    // for ( int i = 0; i < count; i++ ){
    //     for ( int j = 0; j < count; j++ ){
    //         if ( abs(phit(i,j)) < 1e-10)
    //             phit(i,j) = 0.0;
    //     }

    // }

    // cout << "PHI^T*PHI" << endl;
    // cout << phit << endl;
    //     //  phi.col(i) = phi_c.col(i);


    // for (int i = 0; i < count; i++)
    //     phi.col(i) = phi.col(i).normalized();
    
    return phi;

}




MatrixXd SPOD_basis( int Nr, MatrixXd snap, int Nf,  VectorXd &K_pc, VectorXd &lam, MatrixXd &eig_vec, string bc_flag, string filter_flag, double sigma = 1.0){


    int N_snap = snap.cols();

    MatrixXd field(Nr, 1);
    MatrixXd R(N_snap, N_snap);
    MatrixXd R_f(N_snap, N_snap);
    R_f.setZero(N_snap, N_snap);
    MatrixXd phi_c(Nr, N_snap);
    // MatrixXd eig_vec(N_snap, N_snap);
    VectorXd mean(Nr);
    int count;
    VectorXd g(2*Nf+1);

    if ( filter_flag == "BOX"){

        for (int k = 0; k <= 2*Nf; k++ )
            g(k) = 1.0/(2.0*Nf+1.0); 
    }
    else{
        count = 0;
        if ( sigma == 0.0){
            phi_c = POD_basis( Nr, snap, K_pc, lam, eig_vec);
            cout << "sigma = 0 then only POD could be performed" << endl << endl;
            return phi_c;
        }


        for (int k = -Nf; k <= Nf; k++ ){
                g(count) = exp(-k*k/(2.0*sigma*sigma));
            // g(count) = 1.0/sqrt(2.0*M_PI*sigma*sigma)*exp(-k*k/(2.0*sigma*sigma));
            count++;
        }
    }

    mean = snap.rowwise().mean();
    for (int i = 0; i < N_snap; i++)
        snap.col(i) -= mean;

    R = snap.transpose()*snap;

    // cout << "R = \n" << R << endl; 

    //Zero-padded boundary conditions
    if (bc_flag == "ZERO"){

        for (int i = 0; i < N_snap; i++){
            for (int j = 0; j < N_snap; j++){

                count = 0;

                for (int k = -Nf; k <= Nf; k++){
                    if (i + k < 0 || j + k < 0 || i + k >= N_snap || j + k >= N_snap)
                        R_f(i,j) += 0;
                    else
                        R_f(i,j) += g(count)*R(i+k, j+k);
                    
                    count++; 
    
                }
            }
        }
    }

    //Periodic boundary conditions (DA RIVEDERE!!!! non funziona bene)
    if (bc_flag == "PERIODIC"){

        for (int j = 0; j < N_snap; j++){        
            for (int i = j; i < N_snap; i++){

                int index_i = 0;
                int initial_j = abs(j-i);
                int index_j = abs(j-i) ;                
                VectorXd sub_diagonal = VectorXd::Zero(N_snap-initial_j);

                for ( int k = initial_j ; k < (N_snap) ; k++ ) {
                    sub_diagonal(k-initial_j) = R(index_i,index_j);
                    index_i++;
                    index_j++;
                }

                count = 0;
                int l_sub_diag = sub_diagonal.size();

                for (int k = -Nf; k <= Nf; k++) {
                    
                    int w = 0;
                    while ( i+k+w < 0 ) 
                        w += l_sub_diag;

                    while ( i+k+w >= l_sub_diag ) 
                        w -= l_sub_diag;

                    R_f(i,j) += g(count)*sub_diagonal(i+k+w);
                    count++;

                }   
                R_f(j,i) = R_f(i,j);
            //                 cout << "subdiagonal(" << j << ") = \n" << sub_diagonal.transpose() << endl;
            // cin.get();

            }


        }

    }

    // cout << "Filtered correlation matrix" << endl;
    // cout << R_f << endl;

    // cout << "diff R-R_f" << endl;
    // cout << R-R_f << endl;
    


    EigenSolver<MatrixXd> es(R_f); 
    lam = es.eigenvalues().real();
    eig_vec = es.eigenvectors().real();
    eig_sort(lam, eig_vec);

    // cout << "Eigenvalues:\n " << lam << endl << endl;

    double sum = 0;

    for (int i = 0; i < N_snap; i++){
        sum += lam(i)/lam.sum();
        K_pc(i) = sum;
    }

    // cout << "Energy content : \n" << K_pc << endl << endl;
//        cout << "EIG_T*EIG : " << endl;
//    MatrixXd eigcheckort = eig_vec.transpose()*eig_vec;
//     for ( int i = 0; i < N_snap; i++ ){
//         for ( int j = 0; j < N_snap; j++ ){
//             if ( eigcheckort(i,j) < 1e-12 )
//                 eigcheckort(i,j) = 0.0;
//         }

//     } 
//     cout << eigcheckort << endl;
    // MatrixXd data(N_snap, N_snap+1) ;
    // data << lam, eig_vec;

    // for ( int i = 0; i < Ns ; i++){
    //     cout << lam(i) << "\t" << eig_vec.row(i) << endl; 
    // }

    // cout << "eig_vec : \n" << data  << endl;

    double tol = lam(0)*1e-12;
    phi_c = snap*eig_vec;

    count = 0;
    while ( count < lam.size() && abs(lam(count)) > tol)
            count++;

    MatrixXd phi(Nr,count);
    for ( int i = 0 ; i < count ; i++ )
        phi.col(i) = phi_c.col(i)/sqrt(lam(i));

    // MatrixXd phit = phi.transpose()*phi;
    // MatrixXd coef = snap.transpose()*phi*phit.inverse();
// 
// 
    // cout << "The real coefs matrix is : \n \n" << coef << endl << endl << endl;
    // cout << " that is supposed to be equal to eig_vec : \n \n  " << eig_vec << endl << endl << endl;


    // cout << "phi^T*phi" << endl;
    // cout << phi.transpose()*phi << endl;
    
    // phi.col(i) = 1.0/(N_snap*lam(i))*phi_c.col(i);


    // for (int i = 0; i < phi.cols(); i++)
    //     phi.col(i) = phi.col(i).normalized();


    return phi;

}




// MatrixXd Coefs( MatrixXd phi, MatrixXd snap, int Nr, string flag){

//     int N_snap = snap.cols();
//     VectorXd mean(Nr) ;
//     mean = snap.rowwise().mean();
    
//     for (int i = 0; i < N_snap; i++)
//         snap.col(i) -= mean;

//     if ( flag == "POD")
//         return phi.transpose()*snap;
//     else{
//         MatrixXd T = phi.transpose()*phi;
//         MatrixXd pseudo_inv_phi = T.inverse()*phi.transpose();
//         return pseudo_inv_phi*snap;
//     }

// } 


MatrixXcd DMD_basis( int Nr, MatrixXd &snap, VectorXcd &lam, string flag = "STD"){

    int Ns = snap.cols() - 1;
    MatrixXd dmd0_snap = snap.leftCols(Ns);
    MatrixXd dmd1_snap = snap.rightCols(Ns);

    MatrixXd V(Ns, Ns);
    VectorXd lam_POD(Ns), K_pc(Ns); 
    MatrixXd U = POD_basis(Nr, dmd0_snap, K_pc, lam_POD, V, "DMD");
    MatrixXd Sig_inv, Vp;

    if ( Ns > U.cols() ){
        cout << "Ns > U.cols()" << endl;
        cout << " Number of non-zero modes : " << U.cols() << endl;
        Vp = MatrixXd::Zero(Ns, U.cols());
        Vp = V.block(0, 0, Ns, U.cols() );
        Sig_inv = MatrixXd::Zero(U.cols(), U.cols());

        for ( int i = 0; i < U.cols(); i++ )
            Sig_inv(i, i) = 1.0/sqrt(lam_POD(i)); 
    }
    else{
        Vp = MatrixXd::Zero(Ns, Ns);
        Vp = V;
        Sig_inv = MatrixXd::Zero(Ns, Ns);
        for ( int i = 0; i < Ns; i++ )
            Sig_inv(i, i) = 1.0/sqrt(lam_POD(i));

    }

    
    MatrixXd I = MatrixXd::Identity(Ns, Ns);

    cout << "size of V : (" << Vp.rows() << ", " << Vp.cols() << ")" << endl;
    cout << "size of SigInv : (" << Sig_inv.rows() << ", " << Sig_inv.cols() << ")" << endl;
    cout << "size of U : (" << U.rows() << ", " << U.cols() << ")" << endl;

    MatrixXd errM = I - Vp*Sig_inv*U.transpose()*dmd1_snap;
    double err = errM.norm();
// 
    // cout << "Error of DMD in solving the linear problem is : " << err << endl << endl;

    // cout << " U^T : \n " << U.transpose() << endl << endl;
    // cout << " dmd1_snap : \n " << dmd1_snap << endl << endl;
    // cout << " V : \n " << V << endl << endl;
    // cout << " Sig^-1 : \n " << Sig_inv << endl << endl;
    // cout << " Size U: (" << U.rows() << ", " << U.cols() << ") " << endl;
    // cout << " Size dmd1snap: (" << dmd1_snap.rows() << ", " << dmd1_snap.cols() << ") " << endl;
    // cout << " Size Vp: (" << Vp.rows() << ", " << Vp.cols() << ") " << endl;
    // cout << " Size siginv: (" << Sig_inv.rows() << ", " << Sig_inv.cols() << ") " << endl;

    MatrixXd Atilde = U.transpose()*dmd1_snap*Vp*Sig_inv;
    // cout << "Calculated Atilde" << endl;
    // cout << " Atilde : \n " << Atilde << endl << endl;

    EigenSolver<MatrixXd> es(Atilde); 
    lam = es.eigenvalues();
    MatrixXcd eig_vec = es.eigenvectors();

    // cout << "Calculated eigenvalues and eigenvectors" << endl;
    
    MatrixXcd appo = dmd1_snap*Vp*Sig_inv;

    if (flag == "EXACT"){
        MatrixXcd phi(Nr,Vp.cols());

        for (int i = 0; i < Vp.cols(); i++)
            phi.col(i) = 1.0/lam(i)*appo*eig_vec.col(i);
        
        return phi;
    }

    return U*eig_vec;

}

//MatrixXcd Koopman_basis(){   
//}


VectorXcd coefs_dmd ( MatrixXcd phi, VectorXd Init_state){

    //Moore-Penrose pseudo-inverse with Eigen method
    // MatrixXcd pseudo_inv_phi = phi.completeOrthogonalDecomposition().pseudoInverse();
    
    // cout << " Size phiInv : (" << pseudo_inv_phi.rows() << ", " << pseudo_inv_phi.cols() << ")" << endl;
    complex<double> img(0.0,1.0);

    MatrixXcd phiTphi = phi.conjugate().transpose()*phi;
    MatrixXd phiTphi_R = phiTphi.real();
    MatrixXd invM = phiTphi_R.inverse();
    VectorXcd phi_tb = phi.conjugate().transpose()*Init_state;
    VectorXcd cn = phi_tb.transpose()*invM*phi_tb;
    double angle = 0.5*arg(cn(0));
    VectorXcd temp= phi_tb*exp(-img*angle);
    VectorXd br = invM*temp.real();
    VectorXcd b(invM.rows());
//            
    b = br*exp(img*angle);
    cout << "Calculated Coefficients" << endl;

    // VectorXd c = phir.colPivHouseholderQr().solve(Init_state);
 
    // cout << "Inverse of the matrix of modes :\n " << phiTphi.inverse() << endl;  
    // MatrixXcd pseudo_inv_phi = phiTphi.inverse()*phi.transpose();
    // VectorXcd b = phi.colPivHouseholderQr().solve(Init_state);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            


    // int n = pseudo_inv_phi.cols();
    // int m = pseudo_inv_phi.rows();

    // cout << " Pseudo inverse before : \n " << pseudo_inv_phi << endl << endl;

    // for ( int i = 0; i < m; i++){
    //     for ( int j = 0; j < n; j++){
    //         if (abs(pseudo_inv_phi(i,j).real()) < 1e-12 )
    //             pseudo_inv_phi(i,j) = complex<double> (0.0, pseudo_inv_phi(i,j).imag());

    //         if (abs(pseudo_inv_phi(i,j).imag()) < 1e-12 )
    //             pseudo_inv_phi(i,j) = complex<double> (pseudo_inv_phi(i,j).real(), 0.0);

    //     }
    // }
    // cout << " Pseudo inverse after : \n " << pseudo_inv_phi << endl << endl;  

    // cout << " phiTphi : \n" << pseudo_inv_phi << endl << endl; 

    // cout << " Size phi_inv: (" << pseudo_inv_phi.rows() << ", " << pseudo_inv_phi.cols() << ") " << endl;
    // VectorXcd b = pseudo_inv_phi*Init_state; 
    VectorXcd err_ls = Init_state - phi*b;
    cout << " Error of least square complex method in calculating coefficients : " << err_ls.norm() << endl;

    return b;

}


VectorXcd Rec_field_dmd ( int Nr, MatrixXcd phi, VectorXcd lam, VectorXcd coeffs, double t ){

    double dt = 0.001;
    VectorXcd final_state = VectorXcd::Zero(Nr);
    VectorXcd omega(lam.size());

    for ( int i = 0; i < lam.size(); i ++){
        omega(i) = log(lam(i))/dt;

        // if (i == 10){
        // cout << "dmd eigenvalue : " << lam(i) << endl;
        // cout << "Actual eigenvalue : " << omega(i)*dt << endl;
        // }
    }

    for ( int i = 0; i < lam.size(); i++ )
        final_state += exp(omega(i)*t)*coeffs(i)*phi.col(i);
        

    return final_state;

}






VectorXd RBF_Coefs( MatrixXd Coefficients, int Ns, double Dt, double t_star, double t_init = 0.0){


    MatrixXd Coeffs = Coefficients.leftCols(Ns);
    vector<vector<double>> T(Ns,vector<double>(1));
    vector<double> t(1,t_star);
    T[0][0] = t_init;

    for (int i = 1; i < Ns; i++)
        T[i][0] = T[i-1][0]+Dt;

    //Vector of surrogate coefficients
    vector<double> coefs_intrp(Ns);

    // Create surrogates for coefficients
    vector<rbf> surr_coefs;

    RBF_CONSTANTS rbf_const {Dt,0.0};

    for (int i = 0; i < Ns; i++){
        
        vector<double> coefs ;
        for (int j = 0 ; j < Ns ; j++)
            coefs.push_back(Coeffs(i,j));

        // cout << "Coefficients to interpolate : \n" << Coeffs.row(i) << endl;
        // CUBIC, GAUSSIAN, LINEAR, MULTIQUADRATICS


        surr_coefs.push_back( rbf(T, coefs, LINEAR) );
        surr_coefs[i].build();
        surr_coefs[i].evaluate(t,coefs_intrp[i]);

        // cout << "Value of the interpolated coefficient : " << coefs_intrp[i] << endl;
   
    }

    VectorXd coefs_t(Ns);

    for (int i = 0; i < Ns; i++)
        coefs_t(i) = coefs_intrp[i]; 

    //Reconstruction of the field for which RBM has been applied

    return coefs_t;

}




void write_dat(string filename, MatrixXd A){
    
    ofstream data_array;
    data_array.open(filename.c_str());

    for (int i = 0; i < A.rows(); i++){
        for (int j = 0; j < A.cols(); j++)
            data_array << setprecision(12) << scientific << A(i,j) <<  "\t";    

        data_array << endl;
    }


}

#endif

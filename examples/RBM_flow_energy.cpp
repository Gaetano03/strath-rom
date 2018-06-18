#include "stdinclude.hpp"
#include "POD_basis.hpp"
#include "rw_restart.hpp"








int main(int argc, char *argv[] ){


    cout << "\nMain: start" << endl;

    clock_t chrono_begin, chrono_end;
    double comp_time;

    string flag = argv[1];
    int Ns = atoi(argv[2]);
    int sol_freq = 3;
    int m_skip = atoi(argv[3]);
    int read_sol_freq = m_skip*sol_freq;

    string directory = "su2-restart/";
    string filename = directory + "restart_flow.dat";
    int Nr = Ngrid_points( filename, 22 );

    cout << "Number of grid points: " << Nr << endl;

    int Nc = 22;

    vector<int> n_cols = {4,5,6,7};

    VectorXd x(Nr);
    VectorXd y(Nr); 
    VectorXd z(Nr); 
    MatrixXd mat_xyz(Nr, 3);
    MatrixXd field(Nr, 4);
    MatrixXd snap_u(Nr, Ns);
    MatrixXd snap_v(Nr, Ns);
    MatrixXd snap_w(Nr, Ns);
    MatrixXd snap(3*Nr, Ns);
    MatrixXd phi_u(Nr, Ns);
    MatrixXd phi_v(Nr, Ns);
    MatrixXd phi_w(Nr, Ns);
    MatrixXd phi(3*Nr, Ns);
    MatrixXd Coeffs(Ns, Ns);
    MatrixXd S_Coeffs(Ns, Ns);
    VectorXd Rec_field_tstar_POD_u(Nr);
    VectorXd Rec_field_tstar_POD_v(Nr);
    VectorXd Rec_field_tstar_POD_w(Nr);
    VectorXd Rec_field_tstar_SPOD_u(Nr);
    VectorXd Rec_field_tstar_SPOD_v(Nr);
    VectorXd Rec_field_tstar_SPOD_w(Nr);

    string file_in = "restart_flow_";
    string file_temp, file_coef;
    string format = ".dat";
    string file_out = "rbm-restart/POD/PODRec_";

    vector<int> col_xyz = { 1, 2, 3 };
    read_restartDat( filename, col_xyz, 22, mat_xyz );

    x = mat_xyz.col(0);
    y = mat_xyz.col(1);
    z = mat_xyz.col(2);

    int k = 0;
    int n_t_first_snap = 0; 

    double dt = 0.001;
    double Dt = read_sol_freq*dt;
    double t_init = dt*n_t_first_snap + 0.001;

    cout << "Coordinates stored" << endl;

    for( int i = n_t_first_snap; i < n_t_first_snap + read_sol_freq*Ns; i += read_sol_freq ){

        stringstream buffer;
        buffer << setfill('0') << setw(5) << to_string(i);
        file_temp = directory + file_in + buffer.str() + format;
        cout << file_temp << ":\t"; 
 
        read_restartDat( file_temp, n_cols, Nc, field);
        cout << "Complete!" << endl;
        VectorXd rho = field.col(0);
        VectorXd rho_u = field.col(1);
        VectorXd rho_v = field.col(2);
        VectorXd rho_w = field.col(3);
    
        snap_u.col(k) = rho_u.cwiseQuotient(rho);
        snap_v.col(k) = rho_v.cwiseQuotient(rho);
        snap_w.col(k) = rho_w.cwiseQuotient(rho);
        k++;

    }

    snap << snap_u,
            snap_v,
            snap_w;

    VectorXd lam(Ns);
    VectorXd K_pc(Ns);
    VectorXi Nf(51);

    VectorXd mean_u = snap_u.rowwise().mean();
    VectorXd mean_v = snap_v.rowwise().mean();
    VectorXd mean_w = snap_w.rowwise().mean();
    
    
    if ( flag == "POD"){

        chrono_begin = clock();
        phi = POD_basis( 3*Nr, snap, K_pc, lam, Coeffs);
        chrono_end = clock();
        comp_time = ((double)(chrono_end-chrono_begin))/CLOCKS_PER_SEC;
        cout << "Computational time to calculate modes : " << comp_time << endl;

        int Nmod = phi.cols();
        VectorXd t_step(Ns*m_skip);
        phi_u = phi.topRows(Nr);
        phi_v = phi.middleRows(Nr,Nr);
        phi_w = phi.bottomRows(Nr);

        MatrixXd Sig = MatrixXd::Zero(phi.cols(), phi.cols());

        for ( int i = 0; i < phi.cols(); i++ )
            Sig(i, i) = sqrt(lam(i));

        cout << "Number of snapshots : " << Ns << endl;
        cout << "Number of non-zero modes for V = (u, v, w) : " << Nmod << endl;

        vector<string> headers1(7);
        headers1[0] = "\"PointID\"";
        headers1[1] = "\"x\"";
        headers1[2] = "\"y\"";
        headers1[3] = "\"z\"";
        headers1[4] = "\"uPOD\"";
        headers1[5] = "\"vPOD\"";
        headers1[6] = "\"wPOD\"";

        MatrixXd Rec(Nr,3);

        VectorXd err_POD_u(Ns*m_skip);
        VectorXd err_POD_v(Ns*m_skip);
        VectorXd err_POD_w(Ns*m_skip);
        double t = t_init;
        int count = 0;

        for ( int kk = n_t_first_snap; kk < n_t_first_snap + read_sol_freq*Ns; kk += sol_freq){
            t_step(count) = t;

            stringstream dum;
            dum << setfill('0') << setw(5) << to_string(kk);
            file_temp = directory + file_in + dum.str() + format;

            cout << "Reading exact solution  " << file_temp << endl;

            read_restartDat( file_temp, n_cols, Nc, field);

            VectorXd rho = field.col(0);
            VectorXd rho_u = field.col(1);
            VectorXd rho_v = field.col(2);
            VectorXd rho_w = field.col(3);

            VectorXd U_tstar = rho_u.cwiseQuotient(rho);
            VectorXd V_tstar = rho_v.cwiseQuotient(rho);
            VectorXd W_tstar = rho_w.cwiseQuotient(rho);


//             for ( int i = 0; i < n_t.size(); i++){

            chrono_begin = clock();
            VectorXd coefs_interp = RBF_Coefs( Coeffs.transpose(), Nmod, Dt, t, t_init);


//                 for ( int j = 0; j < Nmod; j++){

            Rec_field_tstar_POD_u = phi_u*Sig*coefs_interp.head(Nmod) + mean_u;
            Rec_field_tstar_POD_v = phi_v*Sig*coefs_interp.head(Nmod) + mean_v;
            Rec_field_tstar_POD_w = phi_w*Sig*coefs_interp.head(Nmod) + mean_w;
            chrono_end = clock();
            comp_time = ((double)(chrono_end-chrono_begin))/CLOCKS_PER_SEC;
            cout << "Computational time to calculate reconstruction : " << comp_time << endl;

            VectorXd diff_POD_u = Rec_field_tstar_POD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_POD_u(count) = diff_POD_u.norm()/U_tstar.norm();
            VectorXd diff_POD_v = Rec_field_tstar_POD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_POD_v(count) = diff_POD_v.norm()/V_tstar.norm();
            VectorXd diff_POD_w = Rec_field_tstar_POD_w - W_tstar;
            err_POD_w(count) = diff_POD_w.norm()/W_tstar.norm();
//             }

//         Rec.col(i) = Rec_field_tstar_POD;
//         headers1[i+3] = "\"u"+ to_string(n_t[i])+"\"";

            cout << "Time step  " <<  t << endl;
            cout << "Error for u using all modes : " << err_POD_u(count) << endl;
            cout << "Error for v using all modes : " << err_POD_v(count) << endl;
            cout << "Error for w using all modes : " << err_POD_w(count) << endl;

            stringstream dum7;
            dum7 << setfill('0') << setw(5) << to_string(kk);
            file_temp = file_out + dum7.str() + format;
            cout << "Writing " << file_temp << endl << endl;
            
            Rec << Rec_field_tstar_POD_u, Rec_field_tstar_POD_v, Rec_field_tstar_POD_w; 
            write_restart3D( file_temp, headers1, x, y, z, Rec);
            t += dt*sol_freq;
            count ++;
        }

        MatrixXd err(Ns*m_skip,4);
        err << t_step, err_POD_u, err_POD_v, err_POD_w;
        cout << "Writing err_time.dat" << endl;

        cout << "Main end" << endl;

}else if ( flag == "SPOD" ){
// 


    string flag2 = "truncate";
    int Nfnum = 7;

    MatrixXd cumsum(Ns, Nfnum);
    MatrixXd lambda(Ns, Nfnum);

    double sigma = 4.0;

    VectorXd err_SPOD_u(Ns*Ns);
    VectorXd err_SPOD_v(Ns*Ns);
    VectorXd err_SPOD_w(Ns*Ns);

    int kk = 0;

    Nf(0) = 0;
    
    for ( int i = 1; i <= Nfnum-1; i ++)
        Nf(i) = 6 + Nf(i-1);
    
    // MatrixXd Rec_u(Nr, Nfnum);
    // MatrixXd Rec_v(Nr, Nfnum);
    
    // vector<string> headers1(2*Nfnum+3);
    vector<string> headers1(5);

    headers1[0] = "\"PointID\"";
    headers1[1] = "\"x\"";
    headers1[2] = "\"y\"";


    MatrixXd sn(Nr, 3);
    // cout << "Reading previous solution " << endl;
    // read_restartDat( directory + "su2-snapshots/restart_flow_00000.dat", n_cols, Nc, sn);
    // VectorXd sn_rho = sn.col(0);
    // VectorXd sn_rho_u = sn.col(1);
    // VectorXd sn_rho_v = sn.col(2);
    // VectorXd sn_u = sn_rho_u.cwiseQuotient(sn_rho);
    // VectorXd sn_v = sn_rho_v.cwiseQuotient(sn_rho);
    double tol = 1e-4;

    // vector<MatrixXd> S_phi_Nf_u;
    // vector<MatrixXd> S_phi_Nf_v;
    // vector<MatrixXd> S_phi_Nf_w;
    // vector<MatrixXd> S_Coeffs_Nf;
    VectorXd Nmod(Nfnum);
    VectorXd t_step(Ns*Ns);
    VectorXd N_modes(Ns*Ns);

    vector<string> headers_cum(7);
    vector<string> headers_err(5);

    headers_cum[0] = "\"Nf=0\"";
    headers_cum[1] = "\"Nf=6\"";
    headers_cum[2] = "\"Nf=12\"";
    headers_cum[3] = "\"Nf=18\"";
    headers_cum[4] = "\"Nf=24\"";
    headers_cum[5] = "\"Nf=30\"";
    headers_cum[6] = "\"Nf=36\"";

    headers_err[0] = "\"Time_step\"";
    headers_err[1] = "\"N_modes\"";
    headers_err[2] = "\"err_u\"";
    headers_err[3] = "\"err_v\"";
    headers_err[4] = "\"err_w\"";

    VectorXd diff_SPOD;


    for ( int i = 0; i < Nfnum; i++){

        cout << "-------- SPOD, Nf = " << Nf(i) << " ----------------" << endl;
        chrono_begin = clock();
        MatrixXd S_phi = SPOD_basis( 3*Nr, snap, Nf(i), K_pc, lam, S_Coeffs, "ZERO", "BOX",  sigma );
        chrono_end = clock();
        comp_time = ((double)(chrono_end-chrono_begin))/CLOCKS_PER_SEC;
        cout << endl <<  "Number of snapshots : " << Ns << endl;
        Nmod(i) = S_phi.cols();
        cout << "Number of non-zero modes for V = (u, v, w) : " << Nmod(i) << endl;
        cout << "Computational time to calculate modes : " << comp_time << endl ;

        MatrixXd S_phi_u = S_phi.topRows(Nr);         
        MatrixXd S_phi_v = S_phi.middleRows(Nr,Nr);
        MatrixXd S_phi_w = S_phi.bottomRows(Nr); 

        cumsum.col(i) = K_pc;
        lambda.col(i) = lam;

        double t = t_init;
        int count = 0;

        for ( int kk = n_t_first_snap + 1; kk < n_t_first_snap + read_sol_freq*Ns; kk += read_sol_freq ){

            cout << endl << "                RECONSTRUCTION             " << endl;
            cout << "Time instant : " << t << endl; 


            stringstream dum;
            dum << setfill('0') << setw(5) << to_string(kk);
            file_temp = directory + file_in + dum.str() + format;

            cout << "Reading exact solution  " << file_temp << endl;

            read_restartDat( file_temp, n_cols, Nc, field);

            VectorXd rho = field.col(0);
            VectorXd rho_u = field.col(1);
            VectorXd rho_v = field.col(2);
            VectorXd rho_w = field.col(3);

            VectorXd U_tstar = rho_u.cwiseQuotient(rho);
            VectorXd V_tstar = rho_v.cwiseQuotient(rho);
            VectorXd W_tstar = rho_w.cwiseQuotient(rho);

            VectorXd coefs_interp = RBF_Coefs( S_Coeffs.transpose(), Nmod(i), Dt, t, t_init);

            MatrixXd Sig = MatrixXd::Zero(Nmod(i), Nmod(i));
            for ( int k = 0; k < Nmod(i); k++ )
                Sig(k, k) = sqrt(lambda(k,i));


            for ( int j = 0; j < Nmod(i); j++){
            
                Rec_field_tstar_SPOD_u = S_phi_u.leftCols(j+1)*Sig.topLeftCorner(j+1, j+1)*coefs_interp.head(j+1) + mean_u;
                diff_SPOD = Rec_field_tstar_SPOD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                err_SPOD_u(count) = diff_SPOD.norm()/U_tstar.norm();
                
                Rec_field_tstar_SPOD_v = S_phi_v.leftCols(j+1)*Sig.topLeftCorner(j+1, j+1)*coefs_interp.head(j+1) + mean_v;
                diff_SPOD = Rec_field_tstar_SPOD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                err_SPOD_v(count) = diff_SPOD.norm()/V_tstar.norm();
                
                Rec_field_tstar_SPOD_w = S_phi_w.leftCols(j+1)*Sig.topLeftCorner(j+1, j+1)*coefs_interp.head(j+1) + mean_w;
                diff_SPOD = Rec_field_tstar_SPOD_w - W_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                err_SPOD_w(count) = diff_SPOD.norm()/W_tstar.norm();


                // if ( (j > 2) && (flag2 == "truncate") && ( max(abs(err_SPOD_u(j,i)-err_SPOD_u(j-1,i)),abs(err_SPOD_v(j,i)-err_SPOD_v(j-1,i))) < tol ) &&
                //     ( max(abs(err_SPOD_u(j,i)-err_SPOD_u(j-2,i)),abs(err_SPOD_v(j,i)-err_SPOD_v(j-2,i))) < tol ) ){
                //         Nmod = j+1;
                //         cout << "Number of modes used for Reconstruction: " << Nmod << endl;
                //         cout << "Value of the error for u : " << err_SPOD_u(j,i) << endl;
                //         cout << "Value of the error for v : " << err_SPOD_v(j,i) << endl;

                //         break;
                //     }

                
                t_step(count) = t;
                N_modes(count) = j+1;
                count++;
            }


            if ( Nmod(i) < Ns ){
                for ( int k = 0; k < Ns - Nmod(i); k++){
                    err_SPOD_u(count) = err_SPOD_u(count-1);
                    err_SPOD_v(count) = err_SPOD_v(count-1);
                    err_SPOD_w(count) = err_SPOD_w(count-1);
                    t_step(count) = t;
                    N_modes(count) = Nmod(i) + k + 1;
                    count++;
                }   
            }

            cout << "Error norm of SPOD on u using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_u(count-1) << endl;
            cout << "Error norm of SPOD on v using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_v(count-1) << endl;
            cout << "Error norm of SPOD on w using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_w(count-1) << endl;
        // err_SPOD = MatrixXd::Zero(Nmod, Nfnum);

        // cout << "eigenvalues: \n" << lam << endl<< endl;
        // stringstream buffer;
        // buffer << setfill('0') << setw(5) << to_string(Nf(i));
        // file_temp = directory + "Outpt-ROM/SPOD/SPOD_modes_Nf_" + buffer.str()+format;
        // file_coef = directory + "Outpt-ROM/SPOD/SPOD_coefs_Nf_" + buffer.str()+format;
        // cout << " Writing " << file_temp << " and  coefficients" << "\t";
            t += dt*read_sol_freq;
            cout << endl;
        }
        
        MatrixXd D(Ns*Ns,5);
        D << t_step, N_modes, err_SPOD_u, err_SPOD_v, err_SPOD_w;
        stringstream dum10;
        dum10 << setfill('0') << setw(5) << to_string(Nf(i));
        file_temp = "rbm-restart/SPOD/Error_reconstruction_Vs-Nmodes-t_Nf" + dum10.str() + format;
        
        write_dat( file_temp, D, headers_err);

        cout << endl << endl;
        //MatrixXd Modes(Nr, 2*Ncut);
    }

 
    
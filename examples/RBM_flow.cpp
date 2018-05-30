
#include "stdinclude.hpp"
#include "POD_basis.hpp"
#include "rw_restart.hpp"



int main(int argc, char *argv[] ){


    cout << "\nMain: start" << endl;

    string flag = argv[1];
    int Ns = atoi(argv[2]);
    int sol_freq = 1;
    int m_skip = atoi(argv[3]);
    int read_sol_freq = m_skip*sol_freq;

    string directory = "/home/gaetano/workspace/Simulations/30P30Nunsteady/19deg/";
    string filename = directory + "restart_flow.dat";
    int Nr = Ngrid_points( filename );

    int Nc = 19;

    vector<int> n_cols = {3, 4, 5};

    VectorXd x(Nr);
    VectorXd y(Nr); 
    MatrixXd mat_xy(Nr, 2);
    MatrixXd field(Nr, 3);
    MatrixXd snap_u(Nr, Ns);
    MatrixXd snap_v(Nr, Ns);
    MatrixXd phi_u(Nr, Ns);
    MatrixXd phi_v(Nr, Ns);
    MatrixXd Coeffs_u(Ns, Ns);
    MatrixXd Coeffs_v(Ns, Ns);
    MatrixXd S_Coeffs_u(Ns, Ns);
    MatrixXd S_Coeffs_v(Ns, Ns);
    VectorXd Rec_field_tstar_POD_u(Nr);
    VectorXd Rec_field_tstar_POD_v(Nr);
    VectorXd Rec_field_tstar_SPOD_u(Nr);
    VectorXd Rec_field_tstar_SPOD_v(Nr);

    // string file_in = "su2-snapshots/restart_flow_";
    string file_comp = "restart_flow_"; 
    string file_in = "su2-snapshots/restart_flow_";
    string file_temp, file_coef_u, file_coef_v;
    string format = ".dat";
    string file_out = "Outpt-ROM/SPOD/SPODRec_";

    vector<int> col_xy = { 1, 2 };
    read_restartDat( filename, col_xy, Nc, mat_xy );

    x = mat_xy.col(0);
    y = mat_xy.col(1);

    int k = 0;
    int n_t_first_snap = 0; 
    vector<int> n_t = {1, 3, 5, 7, 9};
    // vector<int> n_t = {2304, 2310, 2316, 2322, 2328};
    double dt = 0.001;
    double Dt = read_sol_freq*dt;
    double t_init = dt*n_t_first_snap;
    vector<double> t_star(n_t.size());



    for ( int i = 0; i < n_t.size(); i++ )
        t_star[i] = dt*n_t[i];

    VectorXd diff_POD, diff_SPOD;

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
    
        snap_u.col(k) = rho_u.cwiseQuotient(rho);
        snap_v.col(k) = rho_v.cwiseQuotient(rho);
        k++;

    }

    VectorXd lam_u(Ns);
    VectorXd lam_v(Ns);
    VectorXd K_pc_u(Ns);
    VectorXd K_pc_v(Ns);

    int Ncut = 10;

    vector<string> headers_mode(2*Ncut+3);

    // headers_mode[0] = "\"PointID\"";
    // headers_mode[1] = "\"x\"";
    // headers_mode[2] = "\"y\"";
    // headers_mode[3] = "\"phi1_u\"";
    // headers_mode[4] = "\"phi2_u\"";
    // headers_mode[5] = "\"phi3_u\"";
    // headers_mode[6] = "\"phi4_u\"";
    // headers_mode[7] = "\"phi1_v\"";
    // headers_mode[8] = "\"phi2_v\"";
    // headers_mode[9] = "\"phi3_v\"";
    // headers_mode[10] = "\"phi4_v\"";


    headers_mode[0] = "\"PointID\"";
    headers_mode[1] = "\"x\"";
    headers_mode[2] = "\"y\"";
    headers_mode[3] = "\"phi1_u\"";
    headers_mode[4] = "\"phi2_u\"";
    headers_mode[5] = "\"phi3_u\"";
    headers_mode[6] = "\"phi4_u\"";
    headers_mode[7] = "\"phi5_u\"";
    headers_mode[8] = "\"phi6_u\"";
    headers_mode[9] = "\"phi7_u\"";
    headers_mode[10] = "\"phi8_u\"";
    headers_mode[11] = "\"phi9_u\"";
    headers_mode[12] = "\"phi10_u\"";
    headers_mode[13] = "\"phi1_v\"";
    headers_mode[14] = "\"phi2_v\"";
    headers_mode[15] = "\"phi3_v\"";
    headers_mode[16] = "\"phi4_v\"";
    headers_mode[17] = "\"phi5_v\"";
    headers_mode[18] = "\"phi6_v\"";
    headers_mode[19] = "\"phi7_v\"";
    headers_mode[20] = "\"phi8_v\"";
    headers_mode[21] = "\"phi9_v\"";
    headers_mode[22] = "\"phi10_v\"";

    
    VectorXi Nf(51);
    MatrixXd Data(Nr,2);
    VectorXd mean_u = snap_u.rowwise().mean();
    VectorXd mean_v = snap_v.rowwise().mean();
    
    
    if ( flag == "POD"){



        phi_u = POD_basis( Nr, snap_u, K_pc_u, lam_u, Coeffs_u);
        phi_v = POD_basis( Nr, snap_v, K_pc_v, lam_v, Coeffs_v);
        int Nmod_u = phi_u.cols();
        int Nmod_v = phi_v.cols();
        VectorXd t_step(Ns*m_skip);

        MatrixXd phi(Nr, 2*Ncut);
        phi << phi_u.leftCols(Ncut), phi_v.leftCols(Ncut) ; 
        write_restart2D(directory + "Outpt-ROM/POD/POD_Modes.dat", headers_mode, x, y, phi);
        write_dat(directory + "Outpt-ROM/POD/POD_coefs_u.dat", Coeffs_u);
        write_dat(directory + "Outpt-ROM/POD/POD_coefs_v.dat", Coeffs_v);

        MatrixXd Sig_u = MatrixXd::Zero(phi_u.cols(), phi_u.cols());
        MatrixXd Sig_v = MatrixXd::Zero(phi_v.cols(), phi_v.cols());

        for ( int i = 0; i < phi_u.cols(); i++ )
            Sig_u(i, i) = sqrt(lam_u(i));
        
        for ( int i = 0; i < phi_v.cols(); i++ )
            Sig_v(i, i) = sqrt(lam_v(i));


        cout << "Number of snapshots : " << Ns << endl;
        cout << "Number of non-zero modes for u : " << Nmod_u << endl;
        cout << "Number of non-zero modes for v : " << Nmod_v << endl << endl;

        vector<string> headers1(5);
//         MatrixXd Rec(Nr, n_t.size());
        headers1[0] = "\"PointID\"";
        headers1[1] = "\"x\"";
        headers1[2] = "\"y\"";
        headers1[3] = "\"uPOD\"";
        headers1[4] = "\"vPOD\"";

        MatrixXd Rec(Nr,2);

        VectorXd err_POD_u(Ns*m_skip);
        VectorXd err_POD_v(Ns*m_skip);
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

            VectorXd U_tstar = rho_u.cwiseQuotient(rho);
            VectorXd V_tstar = rho_v.cwiseQuotient(rho);


//             MatrixXd err_POD(Nmod, n_t.size());


//             for ( int i = 0; i < n_t.size(); i++){


            VectorXd coefs_interp_u = RBF_Coefs( Coeffs_u.transpose(), phi_u, Dt, t, t_init);
            VectorXd coefs_interp_v = RBF_Coefs( Coeffs_v.transpose(), phi_v, Dt, t, t_init);


//                 for ( int j = 0; j < Nmod; j++){

            Rec_field_tstar_POD_u = phi_u*Sig_u*coefs_interp_u.head(Nmod_u) + mean_u;
            Rec_field_tstar_POD_v = phi_v*Sig_v*coefs_interp_v.head(Nmod_v) + mean_v;

            VectorXd diff_POD_u = Rec_field_tstar_POD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_POD_u(count) = diff_POD_u.norm()/U_tstar.norm();
            VectorXd diff_POD_v = Rec_field_tstar_POD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_POD_v(count) = diff_POD_v.norm()/V_tstar.norm();

//             }

//         Rec.col(i) = Rec_field_tstar_POD;
//         headers1[i+3] = "\"u"+ to_string(n_t[i])+"\"";

            cout << "Time step  " <<  t << endl;
            cout << "Error for u using all modes : " << err_POD_u(count) << endl;
            cout << "Error for v using all modes : " << err_POD_v(count) << endl;

//             // MatrixXd fieldn(Nr, 1);

//             // stringstream dum1;
//             // dum1 << setfill('0') << setw(5) << to_string(n_t[i]-5);
//             // file_temp = directory + file_in + dum1.str() + format;
//             // read_restartDat( file_temp, n_cols, Nc, fieldn);
//             // VectorXd u1 = fieldn.col(0);

//             // stringstream dum2;
//             // dum2 << setfill('0') << setw(5) << to_string(n_t[i]+5);
//             // file_temp = directory + file_in + dum2.str() + format;
//             // read_restartDat( file_temp, n_cols, Nc, fieldn);
//             // VectorXd u2 = fieldn.col(0);


//             // VectorXd Rec_interp = u1 + 1.0/Dt*(t_star[i]-t_init)*(u2-u1);
//             // VectorXd err = Rec_interp - field.col(0);
//             // double err_interp = err.norm()/field.col(0).norm();
//             // cout << "Error using linear interpolation : " << err_interp << endl << endl;

// }   256
            stringstream dum7;
            dum7 << setfill('0') << setw(5) << to_string(kk);
            file_temp = directory + file_out + dum7.str() + format;
            cout << "Writing " << file_temp << endl << endl;
            
            Rec << Rec_field_tstar_POD_u, Rec_field_tstar_POD_v; 
            write_restart2D( file_temp, headers1, x, y, Rec);
            t += dt*sol_freq;
            count ++;
        }
        MatrixXd err(Ns*m_skip,3);
        err << t_step, err_POD_u, err_POD_v;
        cout << "Writing err_time.dat" << endl;
        write_dat(directory + "Outpt-ROM/POD/err_time.dat", err);

        cout << "Main end" << endl;

//     write_dat("Error_reconstruction_Nmodes-vs-timestep_v.dat", err_POD);

} else if ( flag == "SPOD" ){
// 
    // Rec_field_tstar_POD = RBF_Coefs( Coeffs.transpose(), phi, Dt, t_star, t_init);
    // Rec_field_tstar_POD += mean;
  
    // VectorXd a(Nmod);
    // for ( int i = 0; i < Nmod; i++)
    //     a(i) = Coeffs(3,i);
    // Rec_field_tstar_POD = phi*a;

    // write_restart("rec_sieber.dat", { "\"PointID\"", "\"x\"", "\"y\"", "\"rec_sieber\""}, x, y, Rec_field_tstar_POD);


    // diff_POD = Rec_field_tstar_POD - field.col(0);
    // err_POD = diff_POD.norm()/field.col(0).norm();
    // cout << endl;
    // cout << "Error norm of POD :" << setprecision(12) << err_POD << endl;
    // cout << endl;

    int Nfnum = 7;
    MatrixXd cumsum_u(Ns, Nfnum);
    MatrixXd lambda_u(Ns, Nfnum);
    MatrixXd cumsum_v(Ns, Nfnum);
    MatrixXd lambda_v(Ns, Nfnum);

    double sigma = 4.0;

    MatrixXd err_SPOD_u(Ns, Nfnum);
    MatrixXd err_SPOD_v(Ns, Nfnum);

    int kk = 0;

    stringstream dum;
    dum << setfill('0') << setw(5) << to_string(n_t[kk]);
    file_temp = directory + file_in + dum.str() + format;

    cout << "Reading exact solution  " << file_temp << endl;

    read_restartDat( file_temp, n_cols, Nc, field);

    VectorXd rho = field.col(0);
    VectorXd rho_u = field.col(1);
    VectorXd rho_v = field.col(2);
    
    VectorXd U_tstar = rho_u.cwiseQuotient(rho);
    VectorXd V_tstar = rho_v.cwiseQuotient(rho);




    Nf(0) = 0;
    
    for ( int i = 1; i <= Nfnum-1; i ++)
        Nf(i) = 10 + Nf(i-1);
    
    MatrixXd Rec_u(Nr, Nfnum);
    MatrixXd Rec_v(Nr, Nfnum);
    
    // vector<string> headers1(2*Nfnum+3);
    vector<string> headers1(5);

    headers1[0] = "\"PointID\"";
    headers1[1] = "\"x\"";
    headers1[2] = "\"y\"";

    MatrixXd dist_modes_u(Ns,Nfnum);
    MatrixXd dist_modes_v(Ns,Nfnum);

    MatrixXd sn(Nr, 3);
    // cout << "Reading previous solution " << endl;
    // read_restartDat( directory + "su2-snapshots/restart_flow_00000.dat", n_cols, Nc, sn);
    // VectorXd sn_rho = sn.col(0);
    // VectorXd sn_rho_u = sn.col(1);
    // VectorXd sn_rho_v = sn.col(2);
    // VectorXd sn_u = sn_rho_u.cwiseQuotient(sn_rho);
    // VectorXd sn_v = sn_rho_v.cwiseQuotient(sn_rho);

    for ( int i = 0; i < Nfnum; i++){


        MatrixXd S_phi_u = SPOD_basis( Nr, snap_u, Nf(i), K_pc_u, lam_u, S_Coeffs_u, "ZERO", "BOX",  sigma );
        MatrixXd S_phi_v = SPOD_basis( Nr, snap_v, Nf(i), K_pc_v, lam_v, S_Coeffs_v, "ZERO", "BOX",  sigma );



        int Nmod_u = S_phi_u.cols();
        int Nmod_v = S_phi_v.cols();

        cumsum_u.col(i) = K_pc_u;
        lambda_u.col(i) = lam_u;
        cumsum_v.col(i) = K_pc_v;
        lambda_v.col(i) = lam_v;

        // err_SPOD = MatrixXd::Zero(Nmod, Nfnum);

        // cout << "eigenvalues: \n" << lam << endl<< endl;

        cout << endl <<  "Number of snapshots : " << Ns << endl;
        cout << "Number of non-zero modes for u : " << Nmod_u << endl;
        cout << "Number of non-zero modes for v : " << Nmod_v << endl;

        stringstream buffer;
        buffer << setfill('0') << setw(5) << to_string(Nf(i));
        file_temp = directory + "Outpt-ROM/SPOD/SPOD_modes_Nf_" + buffer.str()+format;
        file_coef_u = directory + "Outpt-ROM/SPOD/SPOD_coefsu_Nf_" + buffer.str()+format;
        file_coef_v = directory + "Outpt-ROM/SPOD/SPOD_coefsv_Nf_" + buffer.str()+format;
        cout << " Writing " << file_temp << " and  coefficients" << "\t";

        MatrixXd Modes(Nr, 2*Ncut);
        Modes << S_phi_u.leftCols(Ncut), S_phi_v.leftCols(Ncut); 

        write_restart2D(file_temp, headers_mode, x, y, Modes);
        // write_dat(file_coef_u, S_Coeffs_u);
        // write_dat(file_coef_v, S_Coeffs_v);

        cout << "Complete" << endl;

        // S_Coeffs = Coefs( S_phi, snap, Nr, "SPOD");
        VectorXd coefs_interp_u = RBF_Coefs( S_Coeffs_u.transpose(), S_phi_u, Dt, t_star[kk], t_init);
        VectorXd coefs_interp_v = RBF_Coefs( S_Coeffs_v.transpose(), S_phi_v, Dt, t_star[kk], t_init);

        // diff_SPOD = Rec_field_tstar_SPOD - field.col(0);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        // err_SPOD(i) = diff_SPOD.norm()/field.col(0).norm();
        // cout << "Error norm of SPOD for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD(i) << endl;
        cout << endl;

        for ( int j = 0; j < Nmod_u; j++){

            Rec_field_tstar_SPOD_u = S_phi_u.leftCols(j+1)*coefs_interp_u.head(j+1) + mean_u;
            diff_SPOD = Rec_field_tstar_SPOD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_SPOD_u(j,i) = diff_SPOD.norm()/U_tstar.norm();


            // VectorXd dist = sn_u - S_phi_u.col(j);
            // dist_modes_u(j,i) = dist.norm();

            // cout << "error using " << j << "modes : " << err_SPOD(j,i) << endl;
        }
        
        if ( Nmod_u < Ns){
            for ( int k = 0; k < Ns - Nmod_u; k++)
                err_SPOD_u(Nmod_u+k,i) = err_SPOD_u(Nmod_u-1+k,i);   
        }


        // Rec_u.col(i) = Rec_field_tstar_SPOD_u;
        // headers1[i+3] = "\"u"+ to_string(Nf(i))+"\"";
        
        for ( int j = 0; j < Nmod_v; j++){

            Rec_field_tstar_SPOD_v = S_phi_v.leftCols(j+1)*coefs_interp_v.head(j+1) + mean_v;
            diff_SPOD = Rec_field_tstar_SPOD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_SPOD_v(j,i) = diff_SPOD.norm()/V_tstar.norm();
 
            // VectorXd dist = sn_v - S_phi_v.col(j);
            // dist_modes_v(j,i) = dist.norm();
            // cout << "error using " << j << "modes : " << err_SPOD(j,i) << endl;
        }

        if ( Nmod_v < Ns){
            for ( int k = 0; k < Ns - Nmod_v; k++)
                err_SPOD_v(Nmod_v+k,i) = err_SPOD_v(Nmod_v-1+k,i);   
        }        
        
        MatrixXd Rec(Nr, 2);
        Rec << Rec_field_tstar_SPOD_u, Rec_field_tstar_SPOD_v;

        // Rec_v.col(i) = Rec_field_tstar_SPOD_v;
        // headers1[Nfnum+i+3] = "\"v"+ to_string(Nf(i))+"\"";

        cout << "Error norm of SPOD on u using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_u(Nmod_u-1,i) << endl;
        cout << "Error norm of SPOD on v using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_v(Nmod_v-1,i) << endl;

        // cout << "cycle ended" << endl;
    //     if (i == 0){    
    //         MatrixXd Rec(Nr, 2);
    //         Rec << Rec_field_tstar_SPOD_u, Rec_field_tstar_SPOD_v;  
        headers1[3] = "\"uSPOD\"";   
        headers1[4] = "\"vSPOD\"";  
        stringstream buffer1;
        buffer1 << setfill('0') << setw(5) << to_string(Nf(i));
        file_temp = directory + "Outpt-ROM/SPOD/SPODRec_" + to_string(n_t[kk]) + "_" + buffer1.str() + format; 
        write_restart2D( file_temp, headers1, x, y, Rec );
    //     }
    }



    int index_iu, index_ju, index_iv, index_jv; 

    double Err_SPOD_u_min = err_SPOD_u.minCoeff(&index_iu, &index_ju);

    cout << " Minimum error in U reconstruction for the desired time instant is reached for \n ";
    cout << " Nf = " << Nf(index_ju) << endl;
    cout << " Nmodes = " << index_iu + 1 << " out of " << Ns << endl;


    double Err_SPOD_v_min = err_SPOD_v.minCoeff(&index_iv, &index_jv);

    cout << " Minimum error in V reconstruction for the desired time instant is reached for \n ";
    cout << " Nf = " << Nf(index_jv) << endl;
    cout << " Nmodes = " << index_iv + 1 << " out of " << Ns << endl;

    MatrixXd S_phi_u = SPOD_basis( Nr, snap_u, Nf(index_ju), K_pc_u, lam_u, S_Coeffs_u, "ZERO", "BOX",  sigma );
    MatrixXd S_phi_v = SPOD_basis( Nr, snap_v, Nf(index_jv), K_pc_v, lam_v, S_Coeffs_v, "ZERO", "BOX",  sigma );
    VectorXd coefs_interp_u = RBF_Coefs( S_Coeffs_u.transpose(), S_phi_u, Dt, t_star[kk], t_init);
    VectorXd coefs_interp_v = RBF_Coefs( S_Coeffs_v.transpose(), S_phi_v, Dt, t_star[kk], t_init);
    Rec_field_tstar_SPOD_u = S_phi_u.leftCols(index_iu)*coefs_interp_u.head(index_iu) + mean_u;
    Rec_field_tstar_SPOD_v = S_phi_v.leftCols(index_iv)*coefs_interp_u.head(index_iv) + mean_v;
    
    MatrixXd Rec(Nr, 2);
    Rec << Rec_field_tstar_SPOD_u, Rec_field_tstar_SPOD_v;
    headers1[3] = "\"u"+ to_string(Nf(index_ju))+"\"";
    headers1[4] = "\"v"+ to_string(Nf(index_jv))+"\"";

    write_restart2D( directory + "Outpt-ROM/Reconstructed_fields_SPODminerr.dat", headers1, x, y, Rec);

//     // VectorXd Nfd(51);
//     // double count = 0.0;
//     // for ( int i = 0; i <= 50; i ++){
//     //     Nfd(i) = count;
//     //     count += 2.0;
//     // }


//     // Data << Nfd, err_SPOD;
    write_dat( directory + "Outpt-ROM/SPOD/Cumulative_sum_eigenvalues-vs-Nf_u.dat", cumsum_u);
    write_dat( directory + "Outpt-ROM/SPOD/Error_reconstruction_Nmodes-vs-Nf_u.dat", err_SPOD_u);

    write_dat( directory + "Outpt-ROM/SPOD/Cumulative_sum_eigenvalues-vs-Nf_v.dat", cumsum_v);
    write_dat( directory + "Outpt-ROM/SPOD/Error_reconstruction_Nmodes-vs-Nf_v.dat", err_SPOD_v);

    cout << "Main end" << endl;




} else {


// int kk = 0;
int nt = 7;
double t = 0.001*nt;

VectorXd Init_state_u(Nr), Init_state_v(Nr);

// for (int i = 0; i < Ns; i++)
//     snap_u.col(i) -= mean_u;
VectorXcd lam_dmd_u;
MatrixXcd phi_dmd_u = DMD_basis( Nr, snap_u, lam_dmd_u, "EXACT");
VectorXcd lam_dmd_v;
MatrixXcd phi_dmd_v = DMD_basis( Nr, snap_v, lam_dmd_v, "EXACT");

// cout << "Lambda dmd u : " << lam_dmd_u << endl << endl;

stringstream dum1;
dum1 << setfill('0') << setw(5) << to_string(n_t_first_snap);
file_temp = directory + file_in + dum1.str() + format;
// 
cout << "Reading Initial state " << file_temp << endl;
// 
read_restartDat( file_temp, n_cols, Nc, field);
cout << " Complete" << endl << endl;
VectorXd rho = field.col(0);
VectorXd rho_u = field.col(1);
VectorXd rho_v = field.col(2);
// 
Init_state_u = rho_u.cwiseQuotient(rho);
Init_state_v = rho_v.cwiseQuotient(rho);

// cout << " Size phi_u: (" << phi_dmd_u.rows() << ", " << phi_dmd_u.cols() << ") " << endl;
// cout << " Size Init_state_u : " << Init_state_u.size() << endl; 




// cout << "Lambda dmd v : " << lam_dmd_v << endl << endl;
//cout << " Size phi_v: (" << phi_dmd_v.rows() << ", " << phi_dmd_v.cols() << ") " << endl;


// cout << "Fields reconstructed" << endl;

MatrixXd Mode_dmd(Nr, 2*Ncut);
MatrixXd Mode_dmd_im(Nr, 2*Ncut);
MatrixXcd phi_dmd_u_saved = phi_dmd_u.leftCols(Ncut); 
MatrixXcd phi_dmd_v_saved = phi_dmd_v.leftCols(Ncut); 
Mode_dmd << phi_dmd_u_saved.real(), phi_dmd_v_saved.real();
Mode_dmd_im << phi_dmd_u_saved.imag(), phi_dmd_v_saved.imag();
// Mode_dmd << phi_dmd_u.real(), phi_dmd_v.real();
cout << " Mode stored" << endl << endl;

vector<string> headersdmd_mode(2*Ncut+3);
headersdmd_mode[0] = "\"PointID\"";
headersdmd_mode[1] = "\"x\"";
headersdmd_mode[2] = "\"y\"";
headersdmd_mode[3] = "\"phi_u1\"";
headersdmd_mode[4] = "\"phi_u2\"";
headersdmd_mode[5] = "\"phi_u3\"";
headersdmd_mode[6] = "\"phi_u4\"";
headersdmd_mode[7] = "\"phi_u5\"";
headersdmd_mode[8] = "\"phi_u6\"";
headersdmd_mode[9] = "\"phi_u7\"";
headersdmd_mode[10] = "\"phi_u8\"";
headersdmd_mode[11] = "\"phi_u9\"";
headersdmd_mode[12] = "\"phi_u10\"";
headersdmd_mode[13] = "\"phi_u11\"";
headersdmd_mode[14] = "\"phi_u12\"";
headersdmd_mode[15] = "\"phi_u13\"";
headersdmd_mode[16] = "\"phi_u14\"";
headersdmd_mode[17] = "\"phi_u15\"";
headersdmd_mode[18] = "\"phi_u16\"";
headersdmd_mode[19] = "\"phi_u17\"";
headersdmd_mode[20] = "\"phi_u18\"";
headersdmd_mode[21] = "\"phi_u19\"";
headersdmd_mode[22] = "\"phi_u20\"";
headersdmd_mode[23] = "\"phi_u21\"";
headersdmd_mode[24] = "\"phi_u22\"";
headersdmd_mode[25] = "\"phi_u23\"";
headersdmd_mode[26] = "\"phi_u24\"";
headersdmd_mode[27] = "\"phi_u25\"";
headersdmd_mode[28] = "\"phi_u26\"";
headersdmd_mode[29] = "\"phi_u27\"";
headersdmd_mode[30] = "\"phi_u28\"";
headersdmd_mode[31] = "\"phi_u29\"";
headersdmd_mode[32] = "\"phi_u30\"";
headersdmd_mode[33] = "\"phi_u31\"";
headersdmd_mode[34] = "\"phi_u32\"";
headersdmd_mode[35] = "\"phi_u33\"";
headersdmd_mode[36] = "\"phi_u34\"";
headersdmd_mode[37] = "\"phi_u35\"";
headersdmd_mode[38] = "\"phi_v1\"";
headersdmd_mode[39] = "\"phi_v2\"";
headersdmd_mode[40] = "\"phi_v3\"";
headersdmd_mode[41] = "\"phi_v4\"";
headersdmd_mode[42] = "\"phi_v5\"";
headersdmd_mode[43] = "\"phi_v6\"";
headersdmd_mode[44] = "\"phi_v7\"";
headersdmd_mode[45] = "\"phi_v8\"";
headersdmd_mode[46] = "\"phi_v9\"";
headersdmd_mode[47] = "\"phi_v10\"";
headersdmd_mode[48] = "\"phi_v11\"";
headersdmd_mode[49] = "\"phi_v12\"";
headersdmd_mode[50] = "\"phi_v13\"";
headersdmd_mode[51] = "\"phi_v14\"";
headersdmd_mode[52] = "\"phi_v15\"";
headersdmd_mode[53] = "\"phi_v16\"";
headersdmd_mode[54] = "\"phi_v17\"";
headersdmd_mode[55] = "\"phi_v18\"";
headersdmd_mode[56] = "\"phi_v19\"";
headersdmd_mode[57] = "\"phi_v20\"";
headersdmd_mode[58] = "\"phi_v21\"";
headersdmd_mode[59] = "\"phi_v22\"";
headersdmd_mode[60] = "\"phi_v23\"";
headersdmd_mode[61] = "\"phi_v24\"";
headersdmd_mode[62] = "\"phi_v25\"";
headersdmd_mode[63] = "\"phi_v26\"";
headersdmd_mode[64] = "\"phi_v27\"";
headersdmd_mode[65] = "\"phi_v28\"";
headersdmd_mode[66] = "\"phi_v29\"";
headersdmd_mode[67] = "\"phi_v30\"";
headersdmd_mode[68] = "\"phi_v31\"";
headersdmd_mode[69] = "\"phi_v32\"";
headersdmd_mode[70] = "\"phi_v33\"";
headersdmd_mode[71] = "\"phi_v34\"";
headersdmd_mode[72] = "\"phi_v35\"";


// stringstream dum6;
// dum6 << setfill('0') << setw(5) << to_string(6);
// file_temp = directory + file_in + dum6.str() + format;

// cout << " Reading previous state : " << file_temp << endl;
// read_restartDat( file_temp, n_cols, Nc, field);
// cout << " Complete" << endl << endl;

// VectorXd rho_p = field.col(0);
// VectorXd rho_u_p = field.col(1);
// VectorXd rho_v_p = field.col(2);
// // 
// VectorXd Prev_state_u = rho_u_p.cwiseQuotient(rho_p);
// VectorXd Prev_state_v = rho_v_p.cwiseQuotient(rho_p);
write_restart2D(directory + "Outpt-ROM/DMD/ModesDMDreal.dat", headersdmd_mode, x, y, Mode_dmd);
write_restart2D(directory + "Outpt-ROM/DMD/ModesDMDimag.dat", headersdmd_mode, x, y, Mode_dmd_im);
double err = 0.0;
int count = 8;//(Ns-1)*read_sol_freq;
double t_p = 0.0;
double t_in = 0.0;
VectorXcd b_u = coefs_dmd( phi_dmd_u, Init_state_u );
VectorXcd b_v = coefs_dmd( phi_dmd_v, Init_state_v ); 

//while ( err < 0.7 || count > 100 ){

    //count++;
for ( int i = 0; i < Ns*m_skip-1; i ++){

    
    cout << "Time instant : " << t_in << endl;


    stringstream dum5;
    dum5 << setfill('0') << setw(5) << to_string(i+n_t_first_snap);
    file_temp = directory + file_in + dum5.str() + format;

    cout << "Reading exact solution  " << file_temp << endl;

    read_restartDat( file_temp, n_cols, Nc, field);
    rho = field.col(0);
    rho_u = field.col(1);
    rho_v = field.col(2);
    // 
    VectorXd U_tstar = rho_u.cwiseQuotient(rho);
    VectorXd V_tstar = rho_v.cwiseQuotient(rho);

    VectorXcd Rec_dmd_u = Rec_field_dmd ( Nr, phi_dmd_u, lam_dmd_u, b_u, t_in );     
    VectorXcd Rec_dmd_v = Rec_field_dmd ( Nr, phi_dmd_v, lam_dmd_v, b_v, t_in );   

    VectorXd diff_dmd_u = Rec_dmd_u.real() - U_tstar;
    VectorXd diff_dmd_v = Rec_dmd_v.real() - V_tstar;
    double errdmd_u = diff_dmd_u.norm()/U_tstar.norm();
    double errdmd_v = diff_dmd_v.norm()/V_tstar.norm();

    cout << "Error dmd using all modes and the initial state for u: " << errdmd_u << endl;
    cout << "Error dmd using all modes and the initial state for v: " << errdmd_v << endl;

    MatrixXd Rec(Nr, 2);
    Rec << Rec_dmd_u.real(), Rec_dmd_v.real();

    vector<string> headersrec(5);
    headersrec[0] = "\"PointID\"";
    headersrec[1] = "\"x\"";
    headersrec[2] = "\"y\"";
    headersrec[3] = "\"uDMD\"";
    headersrec[4] = "\"vDMD\"";
    stringstream dum4;
    dum4 << setfill('0') << setw(5) << to_string(i);
    file_temp = directory + file_out + dum4.str() + format;
    cout << "Writing " << file_temp << endl << endl;
    write_restart2D(file_temp, headersrec, x, y, Rec);


    t_in += 0.001;

}

    return 0;
// MatrixXd Modes(Nr, 2*Ncut);
// Modes << phi_dmd_u.real().leftCols(Ncut), phi_dmd_v.real().leftCols(Ncut); 

// write_restart2D("ModesDMD.dat", headers_mode, x, y, Modes);


    // cout << " Computing error using the previous instead of the initial state" << endl;



    // cout << " Size phi_u: (" << phi_dmd_u.rows() << ", " << phi_dmd_u.cols() << ") " << endl;
    // cout << " Size Init_state_u : " << Init_state_u.size() << endl; 

    // b_u = coefs_dmd( phi_dmd_u, Prev_state_u );
    // Rec_dmd_u = Rec_field_dmd ( Nr, phi_dmd_u, lam_dmd_u, b_u, t_p );


    // // cout << "Lambda dmd v : " << lam_dmd_v << endl << endl;
    // //cout << " Size phi_v: (" << phi_dmd_v.rows() << ", " << phi_dmd_v.cols() << ") " << endl;

    // b_v = coefs_dmd( phi_dmd_v, Prev_state_v );
    // Rec_dmd_v = Rec_field_dmd ( Nr, phi_dmd_v, lam_dmd_v, b_v, t_p );

    // diff_dmd_u = Rec_dmd_u.real() - U_tstar;
    // diff_dmd_v = Rec_dmd_v.real() - V_tstar;
    // errdmd_u = diff_dmd_u.norm()/U_tstar.norm();
    // errdmd_v = diff_dmd_v.norm()/V_tstar.norm();

    // cout << "Error dmd using all modes and the previous state for u: " << errdmd_u << endl;
    // cout << "Error dmd using all modes and the previous state for v: " << errdmd_v << endl;

    // Rec << Rec_dmd_u.real(), Rec_dmd_v.real();

    // file_temp = directory + "Outpt-ROM/DMDRec_prevState.dat";
    // write_restart2D(file_temp, headersrec, x, y, Rec);

    //err = max(errdmd_u, errdmd_v);
    //t_p +=0.001;
   // }

}


return 0;

}



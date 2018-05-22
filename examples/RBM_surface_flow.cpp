
#include "stdinclude.hpp"
#include "POD_basis.hpp"
#include "rw_restart.hpp"


/* -------ATTENTION-------------
Make a copy of the restart file to read the coordinate in the 
same directory you read the files
*/

int main( int argc, char *argv[] ){



  cout << "\nMain: start" << endl;

    string flag = argv[1];
    int Ns = 70;
    int sol_freq = 5;
    int m_skip = 2;
    int read_sol_freq = m_skip*sol_freq;
    // VectorXi Global_index;

    string directory = "/home/gaetano/workspace/Simulations/30P30Nunsteady/19deg/backupsolution/";
    string filename = directory + "surface_flow.csv";
    int Nr = read_surf_cord_CSV( filename);

    int Nc = 7;

    // vector<int> n_cols = {4, 5, 6};
    vector<int> n_cols = {4, 5, 6};

    VectorXd x(Nr);
    VectorXd y(Nr); 
    MatrixXd mat_xy(Nr, 2);
    MatrixXd field(Nr, 3);
    MatrixXd snap_p(Nr, Ns);
    MatrixXd snap_tx(Nr, Ns);
    MatrixXd snap_ty(Nr, Ns);
    MatrixXd phi_p(Nr, Ns);
    MatrixXd phi_tx(Nr, Ns);
    MatrixXd phi_ty(Nr, Ns);
    MatrixXd Coeffs_p(Ns, Ns);
    MatrixXd Coeffs_tx(Ns, Ns);
    MatrixXd Coeffs_ty(Ns, Ns);
    MatrixXd S_Coeffs_p(Ns, Ns);
    MatrixXd S_Coeffs_tx(Ns, Ns);
    MatrixXd S_Coeffs_ty(Ns, Ns);
    VectorXd Rec_field_tstar_POD_p(Nr);
    VectorXd Rec_field_tstar_POD_tx(Nr);
    VectorXd Rec_field_tstar_POD_ty(Nr);
    VectorXd Rec_field_tstar_SPOD_p(Nr);
    VectorXd Rec_field_tstar_SPOD_tx(Nr);
    VectorXd Rec_field_tstar_SPOD_ty(Nr);

    string file_in = "surface_flow_";
    string file_temp, file_coef;
    string format = ".csv";
    string file_out = "reModel_flow_";

    vector<int> col_xy = { 1, 2 };
    read_CSV( filename, col_xy, Nc, mat_xy );

    x = mat_xy.col(0);
    y = mat_xy.col(1);


    int k = 0;
    int n_t_first_snap = 0; 
    vector<int> n_t = {5, 15, 25, 45, 65, 95, 125, 155, 185, 235, 285, 335, 385, 455, 555, 655};
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
 
        read_CSV( file_temp, n_cols, Nc, field);
        cout << "Complete!" << endl;
        // VectorXd rho_u = field.col(1);
        snap_p.col(k) = field.col(0);
        snap_tx.col(k) = field.col(1);
        snap_ty.col(k) = field.col(2);
        k++;

    }

    VectorXd lam(Ns);
    VectorXd K_pc(Ns);

    int Ncut = 4;

    vector<string> headers(Ncut+3);

    headers[0] = "\"PointID\"";
    headers[1] = "\"x\"";
    headers[2] = "\"y\"";
    headers[3] = "\"phi1\"";
    headers[4] = "\"phi2\"";
    headers[5] = "\"phi3\"";
    headers[6] = "\"phi4\"";
    
    VectorXi Nf(51);
    VectorXd mean_p = snap_p.rowwise().mean();
    VectorXd mean_tx = snap_tx.rowwise().mean();
    VectorXd mean_ty = snap_ty.rowwise().mean();

    if ( flag == "POD"){

        phi_p = POD_basis( Nr, snap_p, K_pc, lam, Coeffs_p);
        phi_tx = POD_basis( Nr, snap_tx, K_pc, lam, Coeffs_tx);
        phi_ty = POD_basis( Nr, snap_ty, K_pc, lam, Coeffs_ty);
        int Nmod_p = phi_p.cols();
        int Nmod_tx = phi_tx.cols();
        int Nmod_ty = phi_ty.cols();

        cout << "Number of snapshots : " << Ns << endl;
        cout << "Number of non-zero modes p: " << Nmod_p << endl;
        cout << "Number of non-zero modes tx: " << Nmod_tx << endl;
        cout << "Number of non-zero modes ty: " << Nmod_ty << endl;


        MatrixXd err_POD_p(Nmod_p, n_t.size());
        MatrixXd err_POD_tx(Nmod_tx, n_t.size());
        MatrixXd err_POD_ty(Nmod_ty, n_t.size());
        MatrixXd Rec_p(Nr, n_t.size());
        MatrixXd Rec_tx(Nr, n_t.size());
        MatrixXd Rec_ty(Nr, n_t.size());

        vector<string> headersp(n_t.size()+3), headerstx(n_t.size()+3), headersty(n_t.size()+3);
        headersp[0] = "\"PointID\"";
        headersp[1] = "\"x_cord\"";
        headersp[2] = "\"y_cord\"";
        headerstx[0] = "\"PointID\"";
        headerstx[1] = "\"x_cord\"";
        headerstx[2] = "\"y_cord\"";
        headersty[0] = "\"PointID\"";
        headersty[1] = "\"x_cord\"";
        headersty[2] = "\"y_cord\"";

        for ( int i = 0; i < n_t.size(); i++){
        
            stringstream dum;
            dum << setfill('0') << setw(5) << to_string(n_t[i]);
            file_temp = directory + file_in + dum.str() + format;

            cout << "Reading exact solution  " << file_temp << endl;

            read_CSV( file_temp, n_cols, Nc, field);

            VectorXd coefs_interp_p = RBF_Coefs( Coeffs_p.transpose(), phi_p, Dt, t_star[i], t_init);
            VectorXd coefs_interp_tx = RBF_Coefs( Coeffs_tx.transpose(), phi_tx, Dt, t_star[i], t_init);
            VectorXd coefs_interp_ty = RBF_Coefs( Coeffs_ty.transpose(), phi_ty, Dt, t_star[i], t_init);


            for ( int j = 0; j < Nmod_p; j++){

                Rec_field_tstar_POD_p = phi_p.leftCols(j+1)*coefs_interp_p.head(j+1) + mean_p;

                diff_POD = Rec_field_tstar_POD_p - field.col(0);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                err_POD_p(j,i) = diff_POD.norm()/field.col(0).norm();

            }

        Rec_p.col(i) = Rec_field_tstar_POD_p;
        headersp[i+3] = "\"p"+ to_string(n_t[i])+"\"";


            for ( int j = 0; j < Nmod_tx; j++){

                Rec_field_tstar_POD_tx = phi_tx.leftCols(j+1)*coefs_interp_tx.head(j+1) + mean_tx;

                diff_POD = Rec_field_tstar_POD_tx - field.col(1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                err_POD_tx(j,i) = diff_POD.norm()/field.col(1).norm();

            }

        Rec_tx.col(i) = Rec_field_tstar_POD_tx;
        headerstx[i+3] = "\"tx"+ to_string(n_t[i])+"\"";

            for ( int j = 0; j < Nmod_ty; j++){

                Rec_field_tstar_POD_ty = phi_ty.leftCols(j+1)*coefs_interp_ty.head(j+1) + mean_ty;

                diff_POD = Rec_field_tstar_POD_ty - field.col(2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                err_POD_ty(j,i) = diff_POD.norm()/field.col(2).norm();

            }

        Rec_ty.col(i) = Rec_field_tstar_POD_ty;
        headersty[i+3] = "\"ty"+ to_string(n_t[i])+"\"";

            cout << "Time step  " <<  n_t[i] << endl;
            cout << "Error p using all modes : " << err_POD_p(Nmod_p-1,i) << endl;
            cout << "Error tx using all modes : " << err_POD_tx(Nmod_tx-1,i) << endl;
            cout << "Error ty using all modes : " << err_POD_ty(Nmod_ty-1,i) << endl;
            // MatrixXd fieldn(Nr, 1);

            // stringstream dum1;
            // dum1 << setfill('0') << setw(5) << to_string(n_t[i]-5);
            // file_temp = directory + file_in + dum1.str() + format;
            // read_restartDat( file_temp, n_cols, Nc, fieldn);
            // VectorXd u1 = fieldn.col(0);

            // stringstream dum2;
            // dum2 << setfill('0') << setw(5) << to_string(n_t[i]+5);
            // file_temp = directory + file_in + dum2.str() + format;
            // read_restartDat( file_temp, n_cols, Nc, fieldn);
            // VectorXd u2 = fieldn.col(0);


            // VectorXd Rec_interp = u1 + 1.0/Dt*(t_star[i]-t_init)*(u2-u1);
            // VectorXd err = Rec_interp - field.col(0);
            // double err_interp = err.norm()/field.col(0).norm();
            // cout << "Error using linear interpolation : " << err_interp << endl << endl;

}   


        write_restart2D("Reconstructed_fields_p.dat", headersp, x, y, Rec_p);
        write_restart2D("Reconstructed_fields_tx.dat", headerstx, x, y, Rec_tx);
        write_restart2D("Reconstructed_fields_ty.dat", headersty, x, y, Rec_ty);
 //   write_dat("Error_reconstruction_Nmodes-vs-timestep_v.dat", err_POD);

        write_dat("Error_reconstruction_Nmodes-vs-p.dat", err_POD_p);
        write_dat("Error_reconstruction_Nmodes-vs-tx.dat", err_POD_tx);
        write_dat("Error_reconstruction_Nmodes-vs-ty.dat", err_POD_ty);



    return 0;
} else{
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

    int Nfnum = 35;
    MatrixXd cumsum_p(Ns, Nfnum);
    MatrixXd lambda_p(Ns, Nfnum);
    MatrixXd cumsum_tx(Ns, Nfnum);
    MatrixXd lambda_tx(Ns, Nfnum);
    MatrixXd cumsum_ty(Ns, Nfnum);
    MatrixXd lambda_ty(Ns, Nfnum);

    double sigma = 10.0;

    MatrixXd err_SPOD_p(Ns, Nfnum);
    MatrixXd err_SPOD_tx(Ns, Nfnum);
    MatrixXd err_SPOD_ty(Ns, Nfnum);

    VectorXd K_pc_p(Ns);
    VectorXd K_pc_tx(Ns);
    VectorXd K_pc_ty(Ns);

    VectorXd lam_p(Ns);
    VectorXd lam_tx(Ns);
    VectorXd lam_ty(Ns);


    int kk = 3;

    stringstream dum;
    dum << setfill('0') << setw(5) << to_string(n_t[kk]);
    file_temp = directory + file_in + dum.str() + format;

    cout << "Reading exact solution  " << file_temp << endl;

    read_restartDat( file_temp, n_cols, Nc, field);


    Nf(0) = 0;
    
    for ( int i = 1; i <= 50; i ++)
        Nf(i) = 2 + Nf(i-1);
    
    MatrixXd Rec_p(Nr, Nfnum);
    MatrixXd Rec_tx(Nr, Nfnum);
    MatrixXd Rec_ty(Nr, Nfnum);    
    vector<string> headers1(35+3);

    headers1[0] = "\"PointID\"";
    headers1[1] = "\"x\"";
    headers1[2] = "\"y\"";



    for ( int i = 0; i < 35; i++){

        sigma = 2.0*i;

        MatrixXd S_phi_p = SPOD_basis( Nr, snap_p, Nf(i), K_pc_p, lam_p, S_Coeffs_p, "ZERO", "BOX",  sigma );
        MatrixXd S_phi_tx = SPOD_basis( Nr, snap_tx, Nf(i), K_pc_tx, lam_tx, S_Coeffs_tx, "ZERO", "BOX",  sigma );
        MatrixXd S_phi_ty = SPOD_basis( Nr, snap_ty, Nf(i), K_pc_ty, lam_ty, S_Coeffs_ty, "ZERO", "BOX",  sigma );

        int Nmod_p = S_phi_p.cols();
        int Nmod_tx = S_phi_tx.cols();
        int Nmod_ty = S_phi_ty.cols();

        cumsum_p.col(i) = K_pc_p;
        lambda_p.col(i) = lam_p;
        cumsum_tx.col(i) = K_pc_tx;
        lambda_tx.col(i) = lam_tx;
        cumsum_ty.col(i) = K_pc_ty;
        lambda_ty.col(i) = lam_ty;
                        

        // err_SPOD = MatrixXd::Zero(Nmod, Nfnum);

        // cout << "eigenvalues: \n" << lam << endl<< endl;

        cout << endl <<  "Number of snapshots : " << Ns << endl;
        cout << "Number of non-zero modes for p : " << Nmod_p << endl;
        cout << "Number of non-zero modes for tx : " << Nmod_tx << endl;
        cout << "Number of non-zero modes for ty : " << Nmod_ty << endl << endl;

        stringstream buffer;
        buffer << setfill('0') << setw(5) << to_string(Nf(i));
        // file_temp = "SPOD_modes_Nf_"+buffer.str()+format;
        // file_coef = "SPOD_coefs_Nf_"+buffer.str()+format;
        // cout << " Writing " << file_temp << " and " << file_coef << "\t";


// 
        // write_restart(file_temp, headers, x, y, S_phi.leftCols(Ncut));
        // write_dat(file_coef, S_Coeffs);

        // cout << "Complete" << endl;

        // S_Coeffs = Coefs( S_phi, snap, Nr, "SPOD");
        VectorXd coefs_interp_p = RBF_Coefs( S_Coeffs_p.transpose(), S_phi_p, Dt, t_star[kk], t_init);
        VectorXd coefs_interp_tx = RBF_Coefs( S_Coeffs_tx.transpose(), S_phi_tx, Dt, t_star[kk], t_init);
        VectorXd coefs_interp_ty = RBF_Coefs( S_Coeffs_ty.transpose(), S_phi_ty, Dt, t_star[kk], t_init);

        Rec_field_tstar_SPOD_p = S_phi_p*coefs_interp_p + mean_p;
        Rec_field_tstar_SPOD_tx = S_phi_tx*coefs_interp_tx + mean_tx;
        Rec_field_tstar_SPOD_ty = S_phi_ty*coefs_interp_ty + mean_ty;
        // diff_SPOD = Rec_field_tstar_SPOD - field.col(0);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        // err_SPOD(i) = diff_SPOD.norm()/field.col(0).norm();
        // cout << "Error norm of SPOD for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD(i) << endl;
        cout << endl;

        for ( int j = 0; j < Nmod_p; j++){

            Rec_field_tstar_SPOD_p = S_phi_p.leftCols(j+1)*coefs_interp_p.head(j+1) + mean_p;

            diff_SPOD = Rec_field_tstar_SPOD_p - field.col(0);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_SPOD_p(j,i) = diff_SPOD.norm()/field.col(0).norm();
            // cout << "error using " << j << "modes : " << err_SPOD(j,i) << endl;
        }

        for ( int j = 0; j < Nmod_tx; j++){

            Rec_field_tstar_SPOD_tx = S_phi_tx.leftCols(j+1)*coefs_interp_tx.head(j+1) + mean_tx;

            diff_SPOD = Rec_field_tstar_SPOD_tx - field.col(1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_SPOD_tx(j,i) = diff_SPOD.norm()/field.col(1).norm();
            // cout << "error using " << j << "modes : " << err_SPOD(j,i) << endl;
        }

        for ( int j = 0; j < Nmod_ty; j++){

            Rec_field_tstar_SPOD_ty = S_phi_ty.leftCols(j+1)*coefs_interp_ty.head(j+1) + mean_ty;

            diff_SPOD = Rec_field_tstar_SPOD_ty - field.col(2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            err_SPOD_ty(j,i) = diff_SPOD.norm()/field.col(2).norm();
            // cout << "error using " << j << "modes : " << err_SPOD(j,i) << endl;
        }


    // Rec.col(i) = Rec_field_tstar_SPOD;
    // headers1[i+3] = "\"u"+ to_string(Nf(i))+"\"";
    cout << "Error norm of SPOD on p using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_p(Nmod_p-1,i) << endl;
    cout << "Error norm of SPOD on tx using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_tx(Nmod_tx-1,i) << endl;
    cout << "Error norm of SPOD on ty using all modes for Nf = " << Nf(i) << " : " << setprecision(12) << err_SPOD_ty(Nmod_ty-1,i) << endl;

    }

        // write_restart("Reconstructed_fields_SPODt5.dat", headers1, x, y, Rec);
// 
//     // VectorXd Nfd(51);
//     // double count = 0.0;
//     // for ( int i = 0; i <= 50; i ++){
//     //     Nfd(i) = count;
//     //     count += 2.0;
//     // }


//     // Data << Nfd, err_SPOD;
    // write_dat("Cumulative_sum_eigenvalues-vs-Nf.dat", cumsum);
    // write_dat("Error_reconstruction_Nmodes-vs-Nf.dat", err_SPOD);

    cout << "Main end" << endl;

    return 0;
}

return 0;

}



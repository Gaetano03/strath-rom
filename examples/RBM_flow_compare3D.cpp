#include "stdinclude.hpp"
#include "POD_basis.hpp"
#include "rw_restart.hpp"



int main(int argc, char *argv[] ){


    cout << "\nMain: start" << endl;

    clock_t chrono_begin, chrono_end;
    double comp_time;


    int Ns = atoi(argv[1]);
    int sol_freq = 3;
    int m_skip = atoi(argv[2]);
    int read_sol_freq = m_skip*sol_freq;

    int Nc = 22;
    string filename = "su2-snapshots/restart_flow.dat";
    int Nr = Ngrid_points( filename, Nc );
    cout << "Number of grid points " <<  Nr <<  endl;

  

    vector<int> n_cols = { 4, 5, 6, 7};

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
    MatrixXd Sphi_u(Nr, Ns);
    MatrixXd Sphi_v(Nr, Ns);
    MatrixXd Sphi_w(Nr, Ns);
    MatrixXd Sphi(3*Nr, Ns);
    MatrixXd Coeffs(Ns, Ns);
    MatrixXd S_Coeffs(Ns, Ns);
    VectorXd Rec_field_tstar_POD_u(Nr);
    VectorXd Rec_field_tstar_POD_v(Nr);
    VectorXd Rec_field_tstar_POD_w(Nr);
    VectorXd Rec_field_tstar_SPOD_u(Nr);
    VectorXd Rec_field_tstar_SPOD_v(Nr);
    VectorXd Rec_field_tstar_SPOD_w(Nr);
    MatrixXd Rec(Nr,15);

    // string file_in = "su2-snapshots/restart_flow_";
    string file_comp = "restart_flow_"; 
    string file_in = "su2-snapshots/restart_flow_";
    string file_temp, file_coef;
    string format = ".dat";
    string file_out = "OutPt-ROM/ComparePODSPODfields";

    vector<int> col_xyz = { 1, 2, 3 };
    read_restartDat( filename, col_xyz, Nc, mat_xyz );

    x = mat_xyz.col(0);
    y = mat_xyz.col(1);
    z = mat_xyz.col(2);


    int kk;
    int n_t_first_snap = 0; 
    double dt = 0.001;
    double Dt = read_sol_freq*dt;
    double t_init = 0.0;
    double t;
    double En;
    double err_POD_u, err_POD_v, err_POD_w, err_SPOD_u, err_SPOD_v, err_SPOD_w;
    double sigma = 4.0;


    cout << "Input energetic content desired in the reconstruction" << endl;
    cin >> En;
    cout << endl;

    cout << "Input the time instant you want to compare(Snap Number)" << endl;
    cin >> kk;
    cout << endl;

    t = dt*kk;

    VectorXi Nf(4);
    Nf(0) = 9;
    Nf(1) = 18;
    Nf(2) = 27;
    Nf(3) = 36;


    VectorXd diff_POD, diff_SPOD;
    double k = 0;

    for( int i = n_t_first_snap; i < n_t_first_snap + read_sol_freq*Ns; i += read_sol_freq ){

        stringstream buffer;
        buffer << setfill('0') << setw(5) << to_string(i);
        file_temp = file_in + buffer.str() + format;
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


    VectorXd lamPOD(Ns);
    VectorXd K_pcPOD(Ns);
    VectorXd lamSPOD(Ns);
    VectorXd K_pcSPOD(Ns);
    VectorXd mean_u = snap_u.rowwise().mean();
    VectorXd mean_v = snap_v.rowwise().mean();
    VectorXd mean_w = snap_w.rowwise().mean();
    
    stringstream dum;
    dum << setfill('0') << setw(5) << to_string(kk);
    file_temp = file_in + dum.str() + format;
    cout << "Reading exact solution  " << file_temp << endl;
    read_restartDat( file_temp, n_cols, Nc, field);
    VectorXd rho = field.col(0);
    VectorXd rho_u = field.col(1);
    VectorXd rho_v = field.col(2);
    VectorXd rho_w = field.col(3);
    VectorXd U_tstar = rho_u.cwiseQuotient(rho);
    VectorXd V_tstar = rho_v.cwiseQuotient(rho);
    VectorXd W_tstar = rho_w.cwiseQuotient(rho);


    phi = POD_basis( 3*Nr, snap, K_pcPOD, lamPOD, Coeffs);
    int NmodPOD = phi.cols();
    cout << "POD modes Calculated" << endl;
    cout << "Number of non-zero modes POD: " << NmodPOD << endl;

    int count = 0;
    double En_contentPOD = 0.0;
    while ( En_contentPOD < En && count < Ns){
        En_contentPOD = K_pcPOD(count);
        count++;
    }

    int N_recPOD = count; 
    phi_u = phi.topRows(Nr);
    phi_v = phi.middleRows(Nr,Nr);
    phi_w = phi.bottomRows(Nr);

    cout << "Number of modes for the desired energy level in POD : " << N_recPOD << endl;

    MatrixXd Sig = MatrixXd::Zero(N_recPOD,N_recPOD);
    for ( int i = 0; i < N_recPOD; i++ )
        Sig(i, i) = sqrt(lamPOD(i));

    VectorXd coefs_interp = RBF_Coefs( Coeffs.transpose(), Ns, Dt, t, t_init);        
    
    Rec_field_tstar_POD_u = phi_u.leftCols(N_recPOD)*Sig*coefs_interp.head(N_recPOD) + mean_u;
    Rec_field_tstar_POD_v = phi_v.leftCols(N_recPOD)*Sig*coefs_interp.head(N_recPOD) + mean_v;
    Rec_field_tstar_POD_w = phi_w.leftCols(N_recPOD)*Sig*coefs_interp.head(N_recPOD) + mean_w;
    
    Rec.col(0) = Rec_field_tstar_POD_u;
    Rec.col(1) = Rec_field_tstar_POD_v;
    Rec.col(2) = Rec_field_tstar_POD_w;

    VectorXd diff_POD_u = Rec_field_tstar_POD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    err_POD_u = diff_POD_u.norm()/U_tstar.norm();
    VectorXd diff_POD_v = Rec_field_tstar_POD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    err_POD_v = diff_POD_v.norm()/V_tstar.norm();
    VectorXd diff_POD_w = Rec_field_tstar_POD_w - W_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    err_POD_w = diff_POD_w.norm()/W_tstar.norm();

    cout << "Time step  " <<  t << endl;
    cout << "Error for u POD using the desired energy level : " << err_POD_u << endl;
    cout << "Error for v POD using the desired energy level : " << err_POD_v << endl;
    cout << "Error for w POD using the desired energy level : " << err_POD_v << endl << endl;



    for ( int r = 0, jj = 3; r < 4; r++){

        cout << "--------- SPOD  Nf = " << Nf(r) << "--------------" << endl;

        chrono_begin = clock();
        Sphi = SPOD_basis( 3*Nr, snap, Nf(r), K_pcSPOD, lamSPOD, S_Coeffs, "ZERO", "BOX",  sigma );
        chrono_end = clock();
        comp_time = ((double)(chrono_end-chrono_begin))/CLOCKS_PER_SEC;
        cout << "Computational time to calculate modes : " << comp_time << endl;

        int NmodSPOD = Sphi.cols();

        cout << "SPOD modes Calculated" << endl;
        cout << "Number of non-zero modes SPOD: " << NmodSPOD << endl;

        count = 0;
        double En_contentSPOD = 0.0;
        while ( En_contentSPOD < En && count < Ns){
            En_contentSPOD = K_pcSPOD(count);
            count++;
        }

        int N_recSPOD = count; 
        Sphi_u = Sphi.topRows(Nr);
        Sphi_v = Sphi.middleRows(Nr,Nr);
        Sphi_w = Sphi.bottomRows(Nr);    

        cout << "Number of modes for the desired energy level in SPOD : " << N_recSPOD << endl;

        MatrixXd SigS = MatrixXd::Zero(N_recSPOD,N_recSPOD);



        for ( int i = 0; i < N_recSPOD; i++ )
            SigS(i, i) = sqrt(lamSPOD(i));


        VectorXd coefs_interpS = RBF_Coefs( S_Coeffs.transpose(), Ns, Dt, t, t_init);

        Rec_field_tstar_SPOD_u = Sphi_u.leftCols(N_recSPOD)*SigS*coefs_interpS.head(N_recSPOD) + mean_u;
        Rec_field_tstar_SPOD_v = Sphi_v.leftCols(N_recSPOD)*SigS*coefs_interpS.head(N_recSPOD) + mean_v;
        Rec_field_tstar_SPOD_w = Sphi_w.leftCols(N_recSPOD)*SigS*coefs_interpS.head(N_recSPOD) + mean_w;

        VectorXd diff_SPOD_u = Rec_field_tstar_SPOD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        err_SPOD_u = diff_SPOD_u.norm()/U_tstar.norm();
        VectorXd diff_SPOD_v = Rec_field_tstar_SPOD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        err_SPOD_v = diff_SPOD_v.norm()/V_tstar.norm();
        VectorXd diff_SPOD_w = Rec_field_tstar_SPOD_w - W_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        err_SPOD_w = diff_SPOD_w.norm()/W_tstar.norm();

        Rec.col(jj++) = Rec_field_tstar_SPOD_u;
        Rec.col(jj++) = Rec_field_tstar_SPOD_v;
        Rec.col(jj++) = Rec_field_tstar_SPOD_w;

        cout << "Time step  " <<  t << endl;
        cout << "Error for u SPOD using the desired energy level : " << err_SPOD_u << endl;
        cout << "Error for v SPOD using the desired energy level : " << err_POD_v << endl;
        cout << "Error for w SPOD using the desired energy level : " << err_SPOD_v << endl << endl;

}

    cout << "Writing Reconstructed solution" << endl;        

    file_temp = file_out + format;
    cout << file_temp << endl << endl;


    vector<string> headers1(18);
    headers1[0] = "\"PointID\"";
    headers1[1] = "\"x\"";
    headers1[2] = "\"y\"";
    headers1[3] = "\"uPOD\"";
    headers1[4] = "\"vPOD\"";
    headers1[5] = "\"wPOD\"";
    headers1[6] = "\"uSPOD9\"";
    headers1[7] = "\"vSPOD9\"";
    headers1[8] = "\"wSPOD9\"";
    headers1[9] = "\"uSPOD18\"";
    headers1[10] = "\"vSPOD18\"";
    headers1[11] = "\"wSPOD18\"";
    headers1[12] = "\"uSPOD27\"";
    headers1[13] = "\"vSPOD27\"";
    headers1[14] = "\"wSPOD27\"";
    headers1[15] = "\"uSPOD36\"";
    headers1[16] = "\"vSPOD36\"";
    headers1[17] = "\"wSPOD36\"";

 
    write_restart3D( file_temp, headers1, x, y, z, Rec);
    
    cout << "Main end" << endl;


}

#include "stdinclude.hpp"
#include "POD_basis.hpp"
#include "rw_restart.hpp"



int main(int argc, char *argv[] ){


    cout << "\nMain: start" << endl;

    int Ns = atoi(argv[1]);
    int sol_freq = 3;
    int m_skip = atoi(argv[2]);
    int read_sol_freq = m_skip*sol_freq;

    string filename = "su2-snapshots/restart_flow.dat";
    int Nr = Ngrid_points( filename );

    int Nc = 19;

    vector<int> n_cols = {3, 4, 5};

    VectorXd x(Nr);
    VectorXd y(Nr); 
    MatrixXd mat_xy(Nr, 2);
    MatrixXd field(Nr, 3);
    MatrixXd snap_u(Nr, Ns);
    MatrixXd snap_v(Nr, Ns);
    MatrixXd snap(2*Nr, Ns);
    MatrixXd phi_u(Nr, Ns);
    MatrixXd phi_v(Nr, Ns);
    MatrixXd phi(2*Nr, Ns);
    MatrixXd Sphi_u(Nr, Ns);
    MatrixXd Sphi_v(Nr, Ns);
    MatrixXd Sphi(2*Nr, Ns);
    MatrixXd Coeffs(Ns, Ns);
    MatrixXd S_Coeffs(Ns, Ns);
    VectorXd Rec_field_tstar_POD_u(Nr);
    VectorXd Rec_field_tstar_POD_v(Nr);
    VectorXd Rec_field_tstar_SPOD_u(Nr);
    VectorXd Rec_field_tstar_SPOD_v(Nr);
    MatrixXd Rec(Nr,8);

    // string file_in = "su2-snapshots/restart_flow_";
    string file_comp = "restart_flow_"; 
    string file_in = "su2-snapshots/restart_flow_";
    string file_temp, file_coef;
    string format = ".dat";
    string file_out = "Outpt-ROM/ComparePOD22SPODfields";

    vector<int> col_xy = { 1, 2 };
    read_restartDat( filename, col_xy, Nc, mat_xy );

    x = mat_xy.col(0);
    y = mat_xy.col(1);

    int kk;
    int n_t_first_snap = 2301; 
    double dt = 0.0015;
    double Dt = read_sol_freq*dt;
    double t_init = dt*(double)n_t_first_snap;
    double t;
    double En;
    double err_POD_u, err_POD_v, err_SPOD_u, err_SPOD_v;
    double sigma = 4.0;


    cout << "Input energetic content desired in the reconstruction" << endl;
    cin >> En;
    cout << endl;

    cout << "Input the time instant you want to compare(Snap Number)" << endl;
    cin >> kk;
    cout << endl;

    t = dt*kk;

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
    
        snap_u.col(k) = rho_u.cwiseQuotient(rho);
        snap_v.col(k) = rho_v.cwiseQuotient(rho);
        k++;

    }

    stringstream dum;
    dum << setfill('0') << setw(5) << to_string(kk);
    file_temp = file_in + dum.str() + format;
    cout << "Reading exact solution  " << file_temp << endl;
    read_restartDat( file_temp, n_cols, Nc, field);
    VectorXd rho = field.col(0);
    VectorXd rho_u = field.col(1);
    VectorXd rho_v = field.col(2);
    VectorXd U_tstar = rho_u.cwiseQuotient(rho);
    VectorXd V_tstar = rho_v.cwiseQuotient(rho);


    snap << snap_u,
            snap_v;


    VectorXd lamPOD(Ns);
    VectorXd K_pcPOD(Ns);
    VectorXd lamSPOD(Ns);
    VectorXd K_pcSPOD(Ns);
    VectorXd mean_u = snap_u.rowwise().mean();
    VectorXd mean_v = snap_v.rowwise().mean();
    

    phi = POD_basis( 2*Nr, snap, K_pcPOD, lamPOD, Coeffs);

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
    phi_v = phi.bottomRows(Nr);

    cout << "Number of modes for the desired energy level in POD : " << N_recPOD << endl;

    MatrixXd Sig = MatrixXd::Zero(N_recPOD,N_recPOD);
    for ( int i = 0; i < N_recPOD; i++ )
        Sig(i, i) = sqrt(lamPOD(i));


    VectorXd coefs_interp = RBF_Coefs( Coeffs.transpose(), Ns, Dt, t, t_init);        
    
    Rec_field_tstar_POD_u = phi_u.leftCols(N_recPOD)*Sig*coefs_interp.head(N_recPOD) + mean_u;
    Rec_field_tstar_POD_v = phi_v.leftCols(N_recPOD)*Sig*coefs_interp.head(N_recPOD) + mean_v;
    Rec.col(0) = Rec_field_tstar_POD_u;
    Rec.col(1) = Rec_field_tstar_POD_v;

    VectorXd diff_POD_u = Rec_field_tstar_POD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    err_POD_u = diff_POD_u.norm()/U_tstar.norm();
    VectorXd diff_POD_v = Rec_field_tstar_POD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    err_POD_v = diff_POD_v.norm()/V_tstar.norm();

    cout << "Time step  " <<  t << endl;
    cout << "Error for u POD using the desired energy level : " << err_POD_u << endl;
    cout << "Error for v POD using the desired energy level : " << err_POD_v << endl << endl;

    VectorXi Nf(3);
    Nf(0) = 10;
    Nf(1) = 20;
    Nf(2) = 30;
    //Nf(3) = 40;


    int Nfnum = 3;

    for ( int i = 0, jj = 2; i < Nfnum; i++ ){

        cout << "---------SPOD Nf =  " << Nf(i) << "---------------" << endl;
        Sphi = SPOD_basis( 2*Nr, snap, Nf(i), K_pcSPOD, lamSPOD, S_Coeffs, "ZERO", "BOX",  sigma );

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
        Sphi_v = Sphi.bottomRows(Nr);    

        cout << "Number of modes for the desired energy level in SPOD : " << N_recSPOD << endl;

        MatrixXd SigS = MatrixXd::Zero(N_recSPOD,N_recSPOD);



        for ( int nval = 0; nval < N_recSPOD; nval++ )
            SigS(nval, nval) = sqrt(lamSPOD(nval));



        VectorXd coefs_interpS = RBF_Coefs( S_Coeffs.transpose(), Ns, Dt, t, t_init);

        Rec_field_tstar_SPOD_u = Sphi_u.leftCols(N_recSPOD)*SigS*coefs_interpS.head(N_recSPOD) + mean_u;
        Rec_field_tstar_SPOD_v = Sphi_v.leftCols(N_recSPOD)*SigS*coefs_interpS.head(N_recSPOD) + mean_v;

        VectorXd diff_SPOD_u = Rec_field_tstar_SPOD_u - U_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        err_SPOD_u = diff_SPOD_u.norm()/U_tstar.norm();
        VectorXd diff_SPOD_v = Rec_field_tstar_SPOD_v - V_tstar;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        err_SPOD_v = diff_SPOD_v.norm()/V_tstar.norm();

        Rec.col(jj++) = Rec_field_tstar_SPOD_u;
        Rec.col(jj++) = Rec_field_tstar_SPOD_v;
        
        cout << "Time step  " <<  t << endl;
        cout << "Error for u SPOD using the desired energy level : " << err_SPOD_u << endl;
        cout << "Error for v SPOD using the desired energy level : " << err_SPOD_v << endl << endl;

    }


    cout << "Writing Reconstructed solution" << endl;        

    vector<string> headers1(11);
    headers1[0] = "\"PointID\"";
    headers1[1] = "\"x\"";
    headers1[2] = "\"y\"";
    headers1[3] = "\"uPOD\"";
    headers1[4] = "\"vPOD\"";
    headers1[5] = "\"uSPOD10\"";
    headers1[6] = "\"vSPOD10\"";
    headers1[7] = "\"uSPOD20\"";
    headers1[8] = "\"vSPOD20\"";
    headers1[9] = "\"uSPOD30\"";
    headers1[10] = "\"vSPOD30\"";
    //headers1[11] = "\"uSPOD40\"";
    //headers1[12] = "\"vSPOD40\"";

    file_temp = file_out + format;
    cout << file_temp << endl << endl;
     
    write_restart2D( file_temp, headers1, x, y, Rec);
    
    cout << "Main end" << endl;


}

#include "stdinclude.hpp"
#include "Vor_detection.hpp"
#include "rw_restart.hpp"

int main( int argc, char *argv[] ){

    printf("Program name %s\n", argv[0]);
    if( argc == 6 ) {
        printf("Number of snaps: %s, Delta snap: %s, First snap: %s, Prob dimension: %s, directory: %s\n", argv[1], argv[2], argv[3], argv[4], argv[5]);
    }
    else if( argc > 6 ) {
        printf("Too many arguments supplied.\n");
    }
    else {
        printf("Four argument expected.\n");
    }


    cout << "Main start" << endl;

    int Ns = atoi( argv[1] );
    int read_sol_freq = atoi( argv[2] );
    int first_snap = atoi( argv[3] );

    string cas = argv[4];
    string directory = argv[5];
    string filename = directory + "restart_flow.dat";

    string file_in = directory + "gradients/gradients_";
    string file_out = directory + "vortices_criteria/Vortex_Criteria_";
    string file_temp_in, file_temp_out;
    string format_in = ".csv";
    string format_out = ".dat";



    unsigned int Nr;
    int Nc;
    VectorXd x(Nr);
    VectorXd y(Nr);
    VectorXd z(Nr);
    vector<int> col_grads;
    vector<string> headers_vortex_det;

    if ( cas == "3D"){

        headers_vortex_det.push_back("\"PointID\"");
        headers_vortex_det.push_back("\"x\"");
        headers_vortex_det.push_back("\"y\"");
        headers_vortex_det.push_back("\"z\"");
        headers_vortex_det.push_back("\"Q_criterion_inc\"");
        headers_vortex_det.push_back("\"Q_criterion_comp\"");
        headers_vortex_det.push_back("\"Q_criterion_inv\"");
        headers_vortex_det.push_back("\"N_k_inc\"");
        headers_vortex_det.push_back("\"N_k_comp\"");
        headers_vortex_det.push_back("\"Lam_ci\"");
        headers_vortex_det.push_back("\"N_trusdell\"");
        headers_vortex_det.push_back("\"Vort_magnitude\"");
        headers_vortex_det.push_back("\"Rortex\"");

        Nr = Ngrid_points( filename, 22 );
        cout << "Number of grid points: " << Nr << endl;

        Nc = 33;

        vector<int> col_xyz = { 1, 2, 3 };
        col_grads = { 21, 22, 23, 24, 25, 26, 27, 28, 29};
        int N_grads;
        N_grads = col_grads.size();


        MatrixXd mat_xyz(Nr, 3);
        // MatrixXd mat_Vortex_Det(Nr, 6);

        cout << "Reading coordinates" << endl;
        read_restartDat ( filename, col_xyz, 22, mat_xyz );

        x = mat_xyz.col(0);
        y = mat_xyz.col(1);
        z = mat_xyz.col(2);

    }


    if ( cas == "2D"){

        headers_vortex_det.push_back("\"PointID\"");
        headers_vortex_det.push_back("\"x\"");
        headers_vortex_det.push_back("\"y\"");
        headers_vortex_det.push_back("\"Q_criterion_inc\"");
        headers_vortex_det.push_back("\"Q_criterion_comp\"");
        headers_vortex_det.push_back("\"Q_criterion_inv\"");
        headers_vortex_det.push_back("\"N_k_inc\"");
        headers_vortex_det.push_back("\"N_k_comp\"");
        headers_vortex_det.push_back("\"Lam_ci\"");
        headers_vortex_det.push_back("\"N_trusdell\"");
        headers_vortex_det.push_back("\"Vort_magnitude\"");
        headers_vortex_det.push_back("\"Rortex\"");


        Nr = Ngrid_points( filename, 19 );
        cout << "Number of grid points: " << Nr << endl;

        Nc = 31;
        


        vector<int> col_xy = { 1, 2 };
        col_grads = { 19, 20, 22, 23};

        int N_grads;
        N_grads = col_grads.size();

        MatrixXd mat_xy(Nr, 2);
        // MatrixXd mat_Vortex_Det(Nr, 6);


        cout << "Reading coordinates" << endl;
        read_restartDat ( filename, col_xy, 19, mat_xy );

        x = mat_xy.col(0);
        y = mat_xy.col(1);

    }


    for ( int i = first_snap; i < (first_snap + read_sol_freq*Ns); i += read_sol_freq ){

    // int i = 240; 

        // Initializing the string with the input file name (default format "restart_flow_05di.dat")
        stringstream buffer;
        buffer << setfill('0') << setw(5) << to_string(i);
        file_temp_in = file_in + buffer.str() + format_in;
        file_temp_out = file_out + buffer.str() + format_out;
        

        MatrixXd mat_Vortex_Det = Vortex_detection ( Nr, Nc, col_grads, file_temp_in ); 
        cout << file_temp_in << endl;

        cout << "Calculated vortex criteria" << endl;

        if ( cas == "2D"){
            write_restart2D(file_temp_out, headers_vortex_det, x, y, mat_Vortex_Det);
        }
        else{
            write_restart3D(file_temp_out, headers_vortex_det, x, y, z, mat_Vortex_Det);
        }

        cout << file_temp_out << endl;

    }


    cout << "Main end\n " << endl;

    return 0;

}
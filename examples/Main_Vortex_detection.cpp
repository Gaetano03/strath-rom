
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
    string format_in = ".dat";
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

    if ( cas == "3DTec"){

        headers_vortex_det.push_back("\"x\"");
        headers_vortex_det.push_back("\"y\"");
        headers_vortex_det.push_back("\"z\"");
        headers_vortex_det.push_back("\"Vort_magnitude\"");
        headers_vortex_det.push_back("\"Rortex\"");
        headers_vortex_det.push_back("\"Mach\"");


        Nr = Ngrid_points( filename, 22 );
        cout << "Number of grid points: " << Nr << endl;

        Nc = 9;

        vector<int> col_xyz = { 1, 2, 3 };
        col_grads = { 0, 1, 2, 3, 4, 5, 6, 7, 8};
        int N_grads;
        N_grads = col_grads.size();


        MatrixXd mat_xyz(Nr, 3);
        // MatrixXd mat_Vortex_Det(Nr, 6);

        cout << "Reading coordinates" << "\t";
        read_restartDat ( filename, col_xyz, 22, mat_xyz );
        cout << "Done" << endl << endl ;
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


        Nr = Ngrid_points( filename, 5 );
        cout << "Number of grid points: " << Nr << endl;

        //Nc = 31;
        Nc = 17;
        


        vector<int> col_xy = { 1, 2 };
        // col_grads = { 19, 20, 22, 23};
        col_grads = { 5, 6, 8, 9};

        int N_grads;
        N_grads = col_grads.size();

        MatrixXd mat_xy(Nr, 2);
        // MatrixXd mat_Vortex_Det(Nr, 6);


        cout << "Reading coordinates" << endl;
        read_restartDat ( filename, col_xy, 5, mat_xy );

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

        cout << "Calculated vortex criteria" << endl;

        if ( cas == "2D"){
            write_restart2D(file_temp_out, headers_vortex_det, x, y, mat_Vortex_Det);
        }
        else if( cas == "3D"){
            write_restart3D(file_temp_out, headers_vortex_det, x, y, z, mat_Vortex_Det);
        }
        else{

            MatrixXd Data(Nr,mat_Vortex_Det.cols()+1);
            
            stringstream buffer1;
            buffer1 << setfill('0') << setw(5) << to_string(i);
            //string file_res_in = "home/gaetano/data_prandtl/simulations/rbm-trap-wing/su2-restart/restart_flow_" + buffer1.str() + format_in;
            string file_res_in = "restart_flow_" + buffer1.str() + format_in;
            cout << "Reading Mach field from : " << file_res_in << "\t";
            MatrixXd MM(Nr,1);
            vector<int> col_M = {15};
            read_restartDat ( file_res_in, col_M, 22, MM );
            cout << "Done" << endl;
            VectorXd Mach = MM.col(0);
            Data << mat_Vortex_Det, Mach;
            cout << "Writing vortices files: " << file_temp_out << "\t"; 
            write_restart3D(file_temp_out, headers_vortex_det, x, y, z, Data);
            cout << "Complete!"<< endl << endl;
        }

        cout << file_temp_out << endl;

    }


    cout << "Main end\n " << endl;

    return 0;

}
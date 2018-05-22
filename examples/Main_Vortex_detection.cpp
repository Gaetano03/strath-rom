
#include "stdinclude.hpp"
#include "Vor_detection.hpp"
#include "rw_restart.hpp"

int main( int argc, char *argv[] ){

printf("Program name %s\n", argv[0]);
if( argc == 5 ) {
    printf("Number of snaps: %s, Delta snap: %s, First snap: %s, directory: %s\n", argv[1], argv[2], argv[3], argv[4]);
}
else if( argc > 5 ) {
    printf("Too many arguments supplied.\n");
}
else {
    printf("Four argument expected.\n");
}


cout << "Main start" << endl;

string directory = argv[4];
string filename = directory + "restart_flow.dat";

unsigned int Nr = Ngrid_points( filename, 22 );

cout << "Number of grid points: " << Nr << endl;

int Nc = 33;
int Ns = atoi( argv[1] );
int read_sol_freq = atoi( argv[2] );
int first_snap = atoi( argv[3] );



vector<int> col_xyz = { 1, 2, 3 };
vector<int> col_grads = { 21, 22, 23, 24, 25, 26, 27, 28, 29};
vector<string> headers_vortex_det(13);
int N_grads;
N_grads = col_grads.size();
VectorXd x(Nr);
VectorXd y(Nr);
VectorXd z(Nr);

MatrixXd mat_xyz(Nr, 3);
// MatrixXd mat_Vortex_Det(Nr, 6);


string file_in = directory + "gradients/gradients_";
string file_out = directory + "vortices_criteria/Vortex_Criteria_";
string file_temp_in, file_temp_out;
string format_in = ".csv";
string format_out = ".dat";

headers_vortex_det[0] = "\"PointID\"";
headers_vortex_det[1] = "\"x\"";
headers_vortex_det[2] = "\"y\"";
headers_vortex_det[3] = "\"z\"";
headers_vortex_det[4] = "\"Q_criterion_inc\"";
headers_vortex_det[5] = "\"Q_criterion_comp\"";
headers_vortex_det[6] = "\"Q_criterion_inv\"";
headers_vortex_det[7] = "\"N_k_inc\"";
headers_vortex_det[8] = "\"N_k_comp\"";
headers_vortex_det[9] = "\"Lam_ci\"";
headers_vortex_det[10] = "\"N_trusdell\"";
headers_vortex_det[11] = "\"Vort_magnitude\"";
headers_vortex_det[12] = "\"Rortex\"";


clock_t chrono_begin, chrono_end;
double comp_time;

chrono_begin = clock();

cout << "Reading coordinates" << endl;
read_restartDat ( filename, col_xyz, 22, mat_xyz );

x = mat_xyz.col(0);
y = mat_xyz.col(1);
z = mat_xyz.col(2);




for ( int i = first_snap; i <= (first_snap + read_sol_freq*Ns); i += read_sol_freq ){

// int i = 240; 

    // Initializing the string with the input file name (default format "restart_flow_05di.dat")
    stringstream buffer;
    buffer << setfill('0') << setw(5) << to_string(i);
    file_temp_in = file_in + buffer.str() + format_in;
    file_temp_out = file_out + buffer.str() + format_out;
    cout << file_temp_in << endl; 

    MatrixXd mat_Vortex_Det = Vortex_detection ( Nr, Nc, col_grads, file_temp_in ); 

    cout << "Calculated vortex criteria" << endl;
    
    write_restart3D(file_temp_out, headers_vortex_det, x, y, z, mat_Vortex_Det);

    cout << file_temp_out << endl;

}


chrono_end = clock();
comp_time = ((double)(chrono_end - chrono_begin))/CLOCKS_PER_SEC;

cout << " Computational time :" << comp_time << endl;

cout << "Main end\n" << endl;

return 0;

}
#ifndef RW_RESTART_HPP
#define RW_RESTART_HPP

#include "stdinclude.hpp"

//Open and read output files from paraview

void read_CSV(string filename, vector<int> n_col, int Nc, MatrixXd &field){


    if ( n_col.size() != field.cols() )
        cout << "ERROR : number of cols must matches with size cols!" << endl;

    ifstream flow_data ;
    flow_data.open( filename );    
   
    if (!flow_data.is_open()){
        cout << "File: " << filename << "not found" << endl;    
        return;
    }

    string line_flow_data ;
    int n_row = 0;

    while ( getline( flow_data, line_flow_data ) ){

        if ( n_row == 0 ){
            n_row++;
            continue;
        }

        RowVectorXd point(n_col.size());
        istringstream iss(line_flow_data) ;
        string token;
        double val;
        int count = 0, c = 0; 

        while( getline( iss, token, ',') ) {

            val = stod(token);

            if ( count == n_col[c]){

                point(c) = val;
                c++;
            }
            count ++;
        }
    
        field.row(n_row-1) = point; 
        n_row++;

    }

    flow_data.close();

}


int read_surf_cord_CSV(string filename){ // VectorXi &Global_index


    ifstream flow_data ;
    flow_data.open( filename );    
   
    if (!flow_data.is_open()){
        cout << "File: " << filename << "not found" << endl;    
        return 0;
    }

    string line_flow_data ;
    int n_row = 0; //, c = 0; 

    while ( getline( flow_data, line_flow_data ) ){



        // if ( n_row == 0 ){
        //     n_row++;
        //     continue;
        // }

        // istringstream iss(line_flow_data) ;
        // string token;
        // int count = 0; 

        // while( getline( iss, token, ',') ) {

        //     if ( count == 0){

        //         Global_index(c) = stoi(token);
        //         c++;
        //     }
        //     count ++;
        // }
     
        n_row++;

    }

    flow_data.close();
    
    return (n_row-1);

}

//typedef Matrix<long double,Dynamic,Dynamic> MatrixXld

// Return the number of grid points from restart files

unsigned int Ngrid_points( string filename, int Nc = 19 ){

    string head;
    ifstream flow_data;
    flow_data.open(filename.c_str());

    if(!flow_data.is_open()){
        cout << "\nFile: " << filename << " not found" << endl;
        return 0;
    }

    unsigned int n_row = 0;

    while(!flow_data.eof()){

        flow_data >> head;
        if( head == "AOA=" )
            break;

        for(unsigned int i = 1; i < Nc; i++){
                flow_data >> head;            
            }

        n_row++;        
    }

    flow_data.close();

    return (n_row-1);

}






//Open and read restart_flow output files from SU2

void read_restartDat(string filename, vector<int> n_col, int Nc, MatrixXd &field){


    int Nr = field.rows();
    ifstream flow_data;

    if ( n_col.size() != field.cols() )
        cout << "ERROR : number of cols of field must matches with size cols!" << endl;

    flow_data.open(filename.c_str());

    if(!flow_data.is_open()){
        cout << "\nFile: " << filename << " not found" << endl;
        return;
    }

    string head;
    long double monnezza;

    for ( int i = 0; i < Nc; i++ ) 
        flow_data >> head;    

    for( int row = 0; row < Nr; row++ ){
        if ( !flow_data.good() )
            break;


        int count = 0;

        for ( int col = 0; col < Nc; col++ ){
            flow_data >> head;

            if ( col == n_col[count] ){
            
                
                monnezza = stold(head);
                
                if ( (monnezza!=0.0) && ((abs(monnezza) < 1e-300) || (abs(monnezza) > 1e300))){
                    cout << " Valore monnezza : " << std::setprecision(17) << std::scientific << monnezza <<  " alla riga : "<< row << endl;
                    monnezza = 0.0;
                }

                field(row, count) = monnezza;

                count++;
            }


        }
    }

    flow_data.close();
}




// Write file with restart_flow structure to be post-processed with SU2_SOL in order to be opened 
// by Paraview or Tecplot
void write_restart2D(string filename, vector<string> headers, VectorXd x, VectorXd y, MatrixXd Var){    

    unsigned int Nr = x.size();
    int Nc = headers.size();

    ofstream flow_data;
    flow_data.open(filename.c_str());

    // Write row of Headers
    for (int i = 0; i < Nc; i++)
        flow_data << headers[i] << " "; 

    flow_data << endl;

    // Write fields
    for (unsigned int i = 0; i < Nr; i++){

        flow_data << i+1 << " ";
        flow_data << setprecision(12) << scientific <<  x(i)  << " ";
        flow_data << setprecision(12) << scientific << y(i)  << " ";

        for (int j = 0; j < Nc-3; j++){

            flow_data << setprecision(12) << scientific << Var(i,j) <<  " ";            

        }

        flow_data << endl;

    }

    // Close file
    flow_data.close();

}


void write_restart3D(string filename, vector<string> headers, VectorXd x, VectorXd y, VectorXd z, MatrixXd Var){    

    unsigned int Nr = x.size();
    int Nc = headers.size();

    ofstream flow_data;
    flow_data.open(filename.c_str());

    // Write row of Headers
    for (int i = 0; i < Nc; i++)
        flow_data << headers[i] << " "; 

    flow_data << endl;

    // Write fields
    for (unsigned int i = 0; i < Nr; i++){

        flow_data << setprecision(12) << scientific <<  x(i)  << " ";
        flow_data << setprecision(12) << scientific << y(i)  << " ";
        flow_data << setprecision(12) << scientific << z(i)  << " ";

        for (int j = 0; j < Nc-3; j++){

            flow_data << setprecision(12) << scientific << Var(i,j) <<  " ";            

        }

        flow_data << endl;

    }

    // Close file
    flow_data.close();

}



// void write_dat(string filename, MatrixXd A){
    
//     ofstream data_array;
//     data_array.open(filename.c_str());

//     for (int i = 0; i < A.cols(); i++){
//         for (int j = 0; j < A.rows(); j++)
//             data_array << A(i,j) << "\t";

//         data_array << endl;
//     }


// }


void read_gradient_dat ( string filename, vector<int> col_grads, int Nc, MatrixXd &field ){


int Nr = field.rows();
    ifstream flow_data;

    if ( col_grads.size() != field.cols() )
        cout << "ERROR : number of cols of field must matches with size cols!" << endl;

    flow_data.open(filename.c_str());

    if(!flow_data.is_open()){
        cout << "\nFile: " << filename << " not found" << endl;
        return;
    }

    string head;
    long double monnezza;

    for( int row = 0; row < Nr; row++ ){
        if ( !flow_data.good() ){
            cout << "Something very bad happened" << endl;
            break;
        }

        int count = 0;

        for ( int col = 0; col < Nc; col++ ){
            flow_data >> head;

            monnezza = stold(head);
                
            if ( (monnezza!=0.0) && ((abs(monnezza) < 1e-300) || (abs(monnezza) > 1e300))){
                cout << " Valore monnezza : " << std::setprecision(17) << std::scientific << monnezza <<  " alla riga : "<< row << endl;
                monnezza = 0.0;
            }
                
            field(row, col) = monnezza;

            }

        cout << field.row(row) << endl;
        cin.get();

        }

    flow_data.close();


}


#endif
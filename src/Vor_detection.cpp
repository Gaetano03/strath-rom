#include "Vor_detection.hpp"

double getRandom(double min, double max)
{
    double before = rand() % (int)max + (int)min;
    double after = (double)rand() / RAND_MAX;
    double result = before + after;
    if (result < min || result > max) {
        result = getRandom(min, max);
    }
    return result;
}


VectorXd GradVbar(VectorXd r, MatrixXd gradV ){
    
    double rx, ry, rz;
    double ux, uy, uz,
            vx, vy, vz,
            wx, wy, wz;
    VectorXd f(3);

    rx = r(0);
    ry = r(1);
    rz = r(2);

    ux = gradV(0,0);
    uy = gradV(0,1);
    uz = gradV(0,2);
    vx = gradV(1,0);
    vy = gradV(1,1);
    vz = gradV(1,2);
    wx = gradV(2,0);
    wy = gradV(2,1);
    wz = gradV(2,2);

    f(0) = rx*rx + ry*ry + rz*rz - 1;
    
    f(1) = - rx*(rx*wx - (ux*(ry*ry + rz*rz + rz))/(rz + 1.0) + (rx*ry*vx)/(rz + 1.0)) 
            - ry*(rx*wy - (uy*(ry*ry + rz*rz + rz))/(rz + 1.0) + (rx*ry*vy)/(rz + 1.0)) 
            - rz*(rx*wz - (uz*(ry*ry + rz*rz + rz))/(rz + 1.0) + (rx*ry*vz)/(rz + 1.0));

    f(2) = - rx*(ry*wx - (vx*(rx*rx + rz*rz + rz))/(rz + 1.0) + (rx*ry*ux)/(rz + 1.0)) 
            - ry*(ry*wy - (vy*(rx*rx + rz*rz + rz))/(rz + 1.0) + (rx*ry*uy)/(rz + 1.0)) 
            - rz*(ry*wz - (vz*(rx*rx + rz*rz + rz))/(rz + 1.0) + (rx*ry*uz)/(rz + 1.0)); 
           

return f;

}
    



MatrixXd Jacobian( VectorXd r, MatrixXd gradV ){

    double rx, ry, rz;
    double ux, uy, uz,
            vx, vy, vz,
            wx, wy, wz;
    MatrixXd J(3,3);

    rx = r(0);
    ry = r(1);
    rz = r(2);

    ux = gradV(0,0);
    uy = gradV(0,1);
    uz = gradV(0,2);
    vx = gradV(1,0);
    vy = gradV(1,1);
    vz = gradV(1,2);
    wx = gradV(2,0);
    wy = gradV(2,1);
    wz = gradV(2,2);


    if (abs(rz + 1.0) < 1e-10){
        rz += getRandom(0.2, 2.0);
        // cout << "Ill-conditioned Jacobian\n changing the initial guess of rz to : " << rz << endl;
        // cin.get();
    }

    J(0,0) = 2*rx;
    
    J(0,1) = 2*ry;
    
    J(0,2) = 2*rz;
    
    J(1,0) = (ux*(ry*ry + rz*rz + rz))/(rz + 1.0) - ry*(wy + (ry*vy)/(rz + 1.0)) 
                - rz*(wz + (ry*vz)/(rz + 1.0)) - rx*wx - rx*(wx + (ry*vx)/(rz + 1.0))
                - (rx*ry*vx)/(rz + 1.0);
    
    J(1,1) = rx*((2.0*ry*ux)/(rz + 1.0) - (rx*vx)/(rz + 1.0)) + ry*((2.0*ry*uy)/(rz + 1.0) 
                - (rx*vy)/(rz + 1.0)) + rz*((2.0*ry*uz)/(rz + 1.0) - (rx*vz)/(rz + 1.0)) 
                - rx*wy + (uy*(ry*ry + rz*rz + rz))/(rz + 1.0) - (rx*ry*vy)/(rz + 1.0);
             
    
    J(1,2) = rx*((ux*(2.0*rz + 1.0))/(rz + 1.0) - (ux*(ry*ry + rz*rz + rz))/pow(rz + 1.0, 2.0) + (rx*ry*vx)/pow(rz + 1.0, 2.0)) 
                - rx*wz + ry*((uy*(2.0*rz + 1.0))/(rz + 1.0) - (uy*(ry*ry + rz*rz + rz))/pow(rz + 1.0, 2.0) 
                + (rx*ry*vy)/pow(rz + 1.0, 2.0)) + rz*((uz*(2.0*rz + 1.0))/(rz + 1.0) - (uz*(ry*ry + rz*rz + rz))/pow(rz + 1.0, 2.0) 
                + (rx*ry*vz)/pow(rz + 1.0, 2.0)) + (uz*(ry*ry + rz*rz + rz))/(rz + 1.0) - (rx*ry*vz)/(rz + 1.0);
    
    J(2,0) = (vx*(rx*rx + rz*rz + rz))/(rz + 1.0) - ry*((ry*uy)/(rz + 1.0) - (2.0*rx*vy)/(rz + 1.0)) - rz*((ry*uz)/(rz + 1.0) 
                - (2.0*rx*vz)/(rz + 1.0)) - ry*wx - rx*((ry*ux)/(rz + 1.0) - (2.0*rx*vx)/(rz + 1.0)) - (rx*ry*ux)/(rz + 1.0);

    J(2,1) = (vy*(rx*rx + rz*rz + rz))/(rz + 1.0) - ry*(wy + (rx*uy)/(rz + 1.0)) - rz*(wz + (rx*uz)/(rz + 1.0)) 
                - ry*wy - rx*(wx + (rx*ux)/(rz + 1.0)) - (rx*ry*uy)/(rz + 1.0);
 

    J(2,2) = rx*((vx*(2.0*rz + 1.0))/(rz + 1.0) - (vx*(rx*rx + rz*rz + rz))/pow(rz + 1.0, 2.0) + (rx*ry*ux)/pow(rz + 1.0, 2.0)) 
                - ry*wz + ry*((vy*(2.0*rz + 1.0))/(rz + 1.0) - (vy*(rx*rx + rz*rz + rz))/pow(rz + 1.0, 2.0) + (rx*ry*uy)/pow(rz + 1.0, 2.0)) 
                + rz*((vz*(2.0*rz + 1.0))/(rz + 1.0) - (vz*(rx*rx + rz*rz + rz))/pow(rz + 1.0, 2.0) + (rx*ry*uz)/pow(rz + 1.0, 2.0)) 
                + (vz*(rx*rx + rz*rz + rz))/(rz + 1.0) - (rx*ry*uz)/(rz + 1.0);
 
 

    return J;

}




VectorXd fsearchzero( VectorXd r0, MatrixXd gradV, double tol, int i, int max_it = 2000, int nt = 0 ){
  
    double eps = 1.0;
    MatrixXd J(3,3);
    VectorXd Dr(3);
    VectorXd f(3);
    VectorXd r = r0;

    int count = 0;

    while ( eps > tol && count < max_it){

        f = -GradVbar( r, gradV );
        J = Jacobian( r, gradV );
        // if ( count == 0)
            // cout << "Jacobian:\n" << J << endl << endl;
        

        Dr = J.colPivHouseholderQr().solve(f);
        
        r += Dr;
        eps = f.norm();  

        // cout << "r :" << setprecision(8) << r.transpose() << "  f  : "<< setprecision(8) << eps << endl;

        count++;
    }

    if ( eps < tol )
        return r;

    // cout << "Convergence not reached" << endl << endl;
    // return r;

    if ( nt < 1 ){
        cout << "Convergence not reached at grid point : " << i << endl;
        cout << "r : " << r << endl;
        cout << "f : " << eps << endl << endl;
        return r;
    }
// 
    // cout << "Other " << nt -1 << "possibilities left" << endl;
    VectorXd rnew(3);
    rnew(0) = getRandom(-1.0, 1.0);
    rnew(1) = getRandom(-1.0, 1.0);
    rnew(2) = getRandom(-1.0, 1.0);

    // rnew = rnew/rnew.norm();
    // cout << "New starting value : " << rnew.transpose() << endl << endl;
    return fsearchzero( rnew, gradV, tol, i, max_it, nt-1 );
        
    
    }

    // cin.get();






MatrixXd Vortex_detection ( unsigned int Nr, int Nc, vector<int> col_grads, string filename ){

    int N_grads = col_grads.size();
    double P, Q, R;
    double Q_crit_inc, Q_crit_c, Q_crit_inv, N_k_inc, N_k_comp, Lam_ci, N_trusdell, Zita, Rortex;
    MatrixXd mat_gradients(Nr, N_grads);


    int dim_grad;    
    if ( N_grads == 4 )
        dim_grad = 2;
    else
        dim_grad = 3;

    // MatrixXd mat_Vortex_Det(Nr, 9);
    MatrixXd mat_Vortex_Det(Nr, 2);
    MatrixXd gradV(dim_grad, dim_grad);
    // MatrixXd Omega(dim_grad, dim_grad);
    // MatrixXd Strain(dim_grad, dim_grad);
    // MatrixXd Strain_D(dim_grad, dim_grad);

    // MatrixXd Omega2(dim_grad, dim_grad);
    // MatrixXd Strain2(dim_grad, dim_grad);
    // MatrixXd Strain_D2(dim_grad, dim_grad);

    VectorXd Vorticity(3);
    VectorXd V_Rortex(3);

    // double n_Omega2, n_Strain2, n_Strain_D2; 
    // double Delta;
    double alfa, beta;

    // VectorXcd lambda(dim_grad);

    cout << "Reading gradient : " << filename <<  endl;
    read_gradient_dat ( filename, col_grads, Nc, mat_gradients );
    cout << "Done" << endl;

    for ( unsigned int i = 0; i < Nr; i++){

        if ( i % 100000 == 0 )
            cout << i << endl;

        int count = 0;

        for ( int m = 0; m < dim_grad; m++ ){
            for ( int n = 0; n < dim_grad; n++ ){
                gradV(m, n) = mat_gradients(i, count);
                count++;
            }
        }

        P = -gradV.trace();        

        if ( dim_grad == 2 ){
            Q = - gradV(1, 0)*gradV(0, 1) + gradV(0, 0)*gradV(1, 1);
            R = 0;
            Zita = abs(gradV(1,0) - gradV(0,1));
        }
        else{
            // Q = -gradV(2, 0)*gradV(0, 2) - gradV(1, 0)*gradV(0, 1)
            //     -gradV(2, 1)*gradV(1, 2) + gradV(1, 1)*gradV(2, 2)
            //     + gradV(0, 0)*gradV(2, 2) + gradV(0, 0)*gradV(1, 1);    
            // R = -gradV.determinant();    
            Vorticity(0) = gradV(2,1) - gradV(1,2);
            Vorticity(1) = gradV(0,2) - gradV(2,0);
            Vorticity(2) = gradV(1,0) - gradV(0,1);
            Zita = Vorticity.norm(); 
        }

        // Omega = 0.5*(gradV - gradV.transpose());
        // Strain = 0.5*(gradV + gradV.transpose());
        // MatrixXd I = MatrixXd::Identity(dim_grad, dim_grad);
        // Strain_D = Strain + P*I;

        // Omega2 = Omega.transpose()*Omega;
        // Strain2 = Strain.transpose()*Strain;
        // Strain_D2 = Strain_D.transpose()*Strain_D;

        // EigenSolver<MatrixXd> eigOmega2(Omega2);
        // EigenSolver<MatrixXd> eigStrain2(Strain2);
        // EigenSolver<MatrixXd> eigStrain_D2(Strain_D2);

        // VectorXcd lam_Omega2(dim_grad);
        // VectorXcd lam_Strain2(dim_grad);
        // VectorXcd lam_Strain_D2(dim_grad);

        // lam_Omega2 = eigOmega2.eigenvalues();
        // lam_Strain2 = eigStrain2.eigenvalues();
        // lam_Strain_D2 = eigStrain_D2.eigenvalues();

        // n_Omega2 = lam_Omega2.real().maxCoeff();
        // n_Strain2 = lam_Strain2.real().maxCoeff();
        // n_Strain_D2 = lam_Strain_D2.real().maxCoeff();

        // Q_crit_inc = 0.5*(n_Omega2 - n_Strain2);
        // Q_crit_c = 0.5*(n_Omega2 - n_Strain_D2);
        // Q_crit_inv = Q;

        // N_k_inc = sqrt(n_Omega2)/sqrt(n_Strain2);
        // N_k_comp = sqrt(n_Omega2)/sqrt(n_Strain_D2);        
        // N_trusdell = Zita/(2.0*Strain_D2.norm()); 

        // Delta = pow( Q/3.0 - pow(P,2.0)/9.0, 3.0) + pow( P*Q/6.0 - pow(P, 3.0)/27.0 - R/2.0, 2.0);

        // if ( Delta > 0 ){
        //     EigenSolver<MatrixXd> eig(gradV);
        //     lambda = eig.eigenvalues();
        //     Lam_ci = lambda.imag().maxCoeff();
        // }
        // else
        //     Lam_ci = 0;


        MatrixXd gradVbar (dim_grad, dim_grad);

        if ( dim_grad == 2 ){
            gradVbar = gradV;

            // cout << "Gradient in rotational frame : \n" << gradVbar << endl << endl;
        }
        else{   //3D case
            MatrixXd Qrot(3,3);
            VectorXd r0(3);
            r0(0) = -1.0;
            r0(1) = 0.0;
            r0(2) = 0.0;

            double tol = 1e-6;
            int max_it = 100;
            double rx, ry, rz;
            
            V_Rortex = fsearchzero( r0, gradV, tol, i, max_it, 1000);
            rx = V_Rortex(0);
            ry = V_Rortex(1);
            rz = V_Rortex(2); 

            // cout << "V_Rortex: \n" << V_Rortex << endl << endl;

            Qrot(0,0) = ( ry*ry +rz*rz +rz )/( 1.0 + rz );
            Qrot(0,1) = -rx*ry/( 1.0 + rz );
            Qrot(0,2) = -rx;
            Qrot(1,0) = -rx*ry/( 1.0 + rz );
            Qrot(1,1) = ( rx*rx +rz*rz + rz )/( 1.0 + rz);
            Qrot(1,2) = -ry;
            Qrot(2,0) = rx;
            Qrot(2,1) = ry;
            Qrot(2,2) = rz;

            gradVbar = Qrot*gradV*Qrot.transpose();
            // cout << "GradVbar : \n" << gradVbar << endl << endl;

        }

        alfa = 0.5*sqrt(pow( gradVbar(1,1) - gradVbar(0,0), 2.0) + pow( gradVbar(1,0) + gradVbar(0,1), 2.0));
        beta = 0.5*( gradVbar(1,0) - gradVbar(0,1) );

        // cout << "alfa, beta:  (" << alfa << "," << beta << ")" << endl << endl;

        if ( (alfa*alfa - beta*beta) > 0 )
            Rortex = 0;
            else if ( beta > 0 ){
                Rortex = 2.0*abs(beta - alfa); // cout << "Rortex diverso da zero" << endl;
            }
                else{ 
                Rortex = 2.0*abs(beta + alfa); // cout << "Rortex diverso da zero" << endl;
                }

        if (isnan(Rortex)){
            cout << "Rortex is Nan at point  " << i << " for some bad reasons" << endl;
            Rortex = 0.0; 

        }

        // cout << "Rortex: " << Rortex << endl;
        // cin.get();


        // mat_Vortex_Det(i,0) = Q_crit_inc;
        // mat_Vortex_Det(i,1) = Q_crit_c;
        // mat_Vortex_Det(i,2) = Q_crit_inv;
        // mat_Vortex_Det(i,3) = N_k_inc;
        // mat_Vortex_Det(i,4) = N_k_comp;
        // mat_Vortex_Det(i,5) = Lam_ci;
        // mat_Vortex_Det(i,6) = N_trusdell;
        // mat_Vortex_Det(i,7) = Zita;
        // mat_Vortex_Det(i,8) = Rortex;
        mat_Vortex_Det(i,0) = Zita;
        mat_Vortex_Det(i,1) = Rortex;

        // cout << "---------------------------------" << endl << endl;
        // cin.get();

    }

    return mat_Vortex_Det;

}



// MatrixXd Rortex_3D ( unsigned int Nr, int Nc, vector<int> col_grads, string filename ){

//     int N_grads = col_grads.size();
//     double Rortex;
//     MatrixXd mat_gradients(Nr, N_grads);


//     int dim_grad = 3;

//     VectorXd M_Rortex(Nr);
//     MatrixXd gradV(dim_grad, dim_grad);

//     VectorXd V_Rortex(3);

//     double Delta;
//     double alfa, beta;

//     read_CSV ( filename, col_grads, Nc, mat_gradients );

//     for ( unsigned int i = 0; i < Nr; i++){

//         cout << i << endl;

//         int count = 0;

//         for ( int m = 0; m < dim_grad; m++ ){
//             for ( int n = 0; n < dim_grad; n++ ){
//                 gradV(m, n) = mat_gradients(i, count);
//                 count++;
//             }
//         }

//         MatrixXd gradVbar (dim_grad, dim_grad);


//         MatrixXd Qrot(3,3);
//         VectorXd r0(3);
//         r0(0) = 1;
//         r0(1) = 1;
//         r0(2) = 1;

//         double tol = 1e-12;
//         double rx, ry, rz;
            
//         V_Rortex = fsearchzero( r0, gradV, tol );
//         rx = V_Rortex(0);
//         ry = V_Rortex(1);
//         rz = V_Rortex(2); 


//         Qrot(0,0) = ( ry*ry +rz*rz +rz )/( 1.0 + rz );
//         Qrot(0,1) = -rx*ry/( 1.0 + rz );
//         Qrot(0,2) = -rx;
//         Qrot(1,0) = -rx*ry/( 1.0 + rz );
//         Qrot(1,1) = ( rx*rx +rz*rz + rz )/( 1.0 + rz);
//         Qrot(1,2) = -ry;
//         Qrot(2,0) = rx;
//         Qrot(2,1) = ry;
//         Qrot(2,2) = rz;

//         gradVbar = Qrot*gradV*Qrot.transpose();

        

//         alfa = 0.5*sqrt(pow( gradVbar(1,1) - gradVbar(0,0), 2.0) + pow( gradVbar(1,0) + gradVbar(0,1), 2.0));
//         beta = 0.5*( gradVbar(1,0) - gradVbar(0,1) );

//         // cout << "alfa, beta:  (" << alfa << "," << beta << ")" << endl << endl;

//         if ( (alfa*alfa - beta*beta) > 0 )
//             Rortex = 0;
//             else if ( beta > 0 )
//                 Rortex = 2.0*abs(beta - alfa);
//                 else 
//                 Rortex = 2.0*abs(beta + alfa);


//         // cout << "Rortex: " << Rortex << endl;
//         // cin.get();

//         M_Rortex(i) = Rortex;

//     }

//     return mat_Vortex_Det;

// }
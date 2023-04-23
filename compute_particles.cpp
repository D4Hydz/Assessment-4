// Assessment 3: n-body gravitational solver

// To avoid warnings tell the compiler to use a recent standard of C++:
// g++ -std=c++17 vector3d.cpp body.cpp compute_orbits.cpp -o compute_orbits

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include "vector3d.hpp"
#include "body.hpp"

using std::cout, std::endl;

// ######################################
// ####### Function declarations ########
// ######################################
// -----------------------------------------------------------
// ADD YOUR FUNCTION DECLARATIONS HERE

// Compute the new velocity and position of each particle in the system 
void update_vel(std::vector<body> &system, double dt);
void update_pos(std::vector<body> &system, double dt);
// Calculate the accelerations of planets in the system at a timestep.
void update_acc(std::vector<body> &system, double epsilon, double sigma, double mass);
// Calculate the new positions and velocites of planets in the system at a given timestep.
void verlet(std::vector<body> &system, double dt);
// Read input data from file
void read_init(std::string input_file, std::vector<body> &system);
// Read the components of a 3d vector from a line
void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z);
// Save the data to file
void save_data(std::ofstream& savefile, const std::vector<body> &system, double t);


// ######################################
// ################ MAIN ################
// ######################################
int main(int argc, char* argv[])
{
  // Checking if number of arguments is equal to 4:
  if (argc != 6) {
    cout << "ERROR: need 4 arguments - compute_orbits <input_file> <output_file> <dt> <T> <Tsave>" << endl;
    return EXIT_FAILURE;
  }
  // Process command line inputs:
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  double dt = atof(argv[3]); // Time step
  int T = atoi(argv[4]); // Total number of time steps
  int Tsave = atoi(argv[5]); // Number of steps between each save

  std::vector<body> system; // Create an empty vector container for bodies
  read_init(input_file, system); // Read bodies from input file into system
  int N = system.size(); // Number of bodies in system

  cout << "--- Orbital motion simulation ---" << endl;
  cout << " number of bodies N: " << N << endl;
  for(int p=0; p < N; p++)
  {
    cout << "- " << system[p].get_name() << endl; // Display names
  }
  cout << "       time step dt: " << dt << endl;
  cout << "  number of steps T: " << T << endl;
  cout << "   save steps Tsave: " << Tsave << endl;

  std::ofstream savefile (output_file); // Open save file
  if (!savefile.is_open()) {
    cout << "Unable to open file: " << output_file << endl; // Exit if save file has not opened
    return EXIT_FAILURE;
  }
  savefile << std::setprecision(16); // Set the precision of the output to that of a double
  // Write a header for the save file
  savefile << dt << "," << T << "," << Tsave << "," << N << endl;
  for(int p=0; p < (N-1); p++)
  {
    savefile << system[p].get_name() << ",";
  }
  savefile << system[N-1].get_name() << endl;

  // -----------------------------------------------------------
  // ADD YOUR CODE HERE
  
  // Initialise variables
  double epsilon = 
  double total_time = 0.;
  int Tsave_counter = 0;
  
  // Calculate and save data for timestep t = 0. 
  update_acc(system, epsilon, sigma, mass); // Calculate acceleration at t = 0
  save_data(savefile, system, total_time);
    
  // Loop steps until it reaches total number for time steps 'T'.
  for (int t=1; t<=T; t++)
  {
    // Runs Vel_verlet function.
    vel_verlet(system, dt);
    
    // Adds 1 to Tsave step counter.
    Tsave_counter++;

    // Checks if step counter is equal to Tsave.
    if (Tsave_counter == Tsave)
    {
      // Calculates the total time at this timestep.
      total_time = t * dt;

      // Computes energies and momentum and saves data to file.
      save_data(savefile, system, total_time);
      
      // Reset Tsave counter to 0.
      Tsave_counter = 0;   
    }
  }
  // -----------------------------------------------------------
  savefile.close();
  return EXIT_SUCCESS; 
}


// ######################################
// ###### Function implementations ######
// ######################################
// -----------------------------------------------------------
// ADD YOUR FUNCTION IMPLEMENTATIONS HERE
void verlet(std::vector<body> &system, double dt, double epsilon, double sigma, double mass)
{
  // Completes the leap frog integration scheme in order.
  update_vel(system, dt);
  update_pos(system, dt);
  update_acc(system, epsilon, sigma, mass);
  update_vel(system, dt);
} 

void update_pos(std::vector<body> &system, double dt)
{ 
  int N = system.size();
  for (int i = 0; i < N; i++)
  {
    // Initialise variables for this planet.
    vec pos = system[i].get_pos();
    vec velocity = system[i].get_vel();

    pos += (dt*(velocity)); // Calculates new position based on new velocity.
    system[i].set_pos(pos); // Sets new position.
  }
}

void update_vel(std::vector<body> &system, double dt)
{ 
  int N = system.size();
  for (int i = 0; i < N; i++)
  {
    // Initialise variables for this planet.
    vec velocity = system[i].get_vel();
    vec acceleration = system[i].get_acc();
    
    velocity +=(dt / 2 *(acceleration)); // Calculates new velocity over half a timestep.
    system[i].set_vel(velocity); // Sets new velocity to particle.
  }
}

void update_acc(std::vector<body> &system, double epsilon, double sigma, double mass)
{
  int N = system.size();

  for (int i=0; i<N; i++) // Loop through particles.
    {
      // Initialise variables for first particle.
      vec tot_force(0.,0.,0.);
      vec posi = system[i].get_pos();
      
      for (int j=0; j<N; j++) // Loop through other particles.
      {
        if (i != j) // Particles must be different!
        {
          vec posj = system[j].get_pos(); // Initialise variables for second particle.
          
          vec distance = posj - (posi); // Calculates the vector distance between them.
          
          double length = distance.length(); // Calculates the length of the vector.

          // Calculates the scalar force component between the 2 particles
          vec force = distance / (length) * (24 * (epsilon / length) * (2 * pow((sigma / length), 12) - pow((sigma / length), 6)));
          
          tot_force += (force); // Sums the different forces on the particle.
        }
      }
      system[i].set_acc(tot_force / (mass)); // Calculates and sets acceleration variable of planet.
    }
}

// -----------------------------------------------------------
void read_init(std::string input_file, std::vector<body> &system)
{
  std::string line; // Declare a string to store each line
  std::string name; // String to store body name
  double m, x, y, z, vx, vy, vz; // Doubles to store vector components
  int line_cnt = 0; // Line counter

  // Declare and initialise an input file stream object
  std::ifstream data_file(input_file); 

  while (getline(data_file, line)) // Read the file line by line
  {
    line_cnt++;
    std::stringstream data_line(line); // Create a string stream from the line
    switch (line_cnt)
    {
      case 1:
        name = line;
        break;
      case 2:
        m = std::stod(line); // Convert string line into double
        break;            
      case 3:
        read_vector3d(data_line, x, y, z); // Read the 3 components of the vector on the line
        break;
      case 4:
        read_vector3d(data_line, vx, vy, vz); // Read the 3 components of the vector on the line
        break;                   
      }
      if (line_cnt==4) // Data for one body has been extracted
      {
        line_cnt = 0; // Reset line counter
        body b(name,m,vec(x,y,z),vec(vx,vy,vz)); // Package data into body
        system.push_back(b); // Add body to system
      }
    }
    // Close the file
    data_file.close();
}

void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z)
{
  std::string value; // Declare a string to store values in a line
  int val_cnt = 0; // Value counter along the line
  
  while (getline(data_line, value, ','))
  {
    val_cnt++;
    switch (val_cnt)
    {
      case 1:
        x = std::stod(value); 
        break;
      case 2:
        y = std::stod(value); 
        break;
      case 3:
        z = std::stod(value);
        break;
    }
  } 
}

void save_data(std::ofstream& savefile, const std::vector<body> &system, double t)
{
    // Function for saving the simulation data to file.

    vec L; // Total angular momentum
    double E = 0.0; // Total energy
    for(int p = 0; p < system.size(); p++)
    { 
      E += system[p].get_ke() + 0.5*system[p].get_gpe();
      L += system[p].get_L();
    }
    double Lmag = L.length(); // Magnitude of total angular momentum

    // Write a header for this time-step with time, total energy and total mag of L:
    savefile << t << "," << E << "," << Lmag << endl;

    // Loop over the bodies:
    for(int p = 0; p < system.size(); p++)
    { 
      // Output position and velocity for each body:
      savefile << system[p].get_pos() << "," << system[p].get_vel() << endl;
    }
}
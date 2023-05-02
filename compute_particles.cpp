// Assessment 4: N-body Molecular Dynamics Model

// g++ -std=c++17 vector3d.cpp body.cpp compute_particles.cpp -o compute_particles
// ./compute_particles particles.csv test.csv 1e-5 1080 1

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

// Compute the new velocity and position of each particle in the system.
void update_vel(std::vector<body> &system, double dt);
void update_pos(std::vector<body> &system, double dt, double boundry);
// Calculate the accelerations of planets in the system at a timestep.
void update_acc(std::vector<body> &system, double epsilon, double sigma, double mass);
// Calculates the temperature of the system.
double calc_temp(std::vector<body> &system, double mass);
// Read input data from file
void read_init(std::string input_file, std::vector<body> &system);
// Read the components of a 3d vector from a line
void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z);
// Save the data to file
void save_data(std::ofstream& savefile, const std::vector<body> &system);

// ######################################
// ################ MAIN ################
// ######################################
int main(int argc, char* argv[])
{
  // Checking if number of arguments is equal to 4:
  if (argc != 7) {
    cout << "ERROR: need 6 arguments - compute_particles <input_file> <output_file> <dt> <T> <Tsave> <boundry>" << endl;
    return EXIT_FAILURE;
  }
  // Process command line inputs:
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  double dt = atof(argv[3]); // Time step
  int T = atoi(argv[4]); // Total number of time steps
  int Tsave = atoi(argv[5]); // Number of steps between each save
  double boundry = atoi(argv[6]); // Size of boundry limit in angstroms (Å)

  std::vector<body> system; // Create an empty vector container for bodies
  read_init(input_file, system); // Read bodies from input file into system
  int N = system.size(); // Number of bodies in system

  cout << "--- Molecular Dynamics simulation ---" << endl;
  cout << " number of bodies N: " << N << endl;
  cout << "       time step dt: " << dt << endl;
  cout << "  number of steps T: " << T << endl;
  cout << "   save steps Tsave: " << Tsave << endl;
  cout << "half boundry length: " << boundry << endl;
  
  // for(int p=0; p < N; p++) // Display initial positions and velocities.
  // {
  //   cout << "Particle " << p << " position " << system[p].get_pos() << endl;
  //   cout << "Particle " << p << " velocity " << system[p].get_vel() << endl;
  // }

  std::ofstream savefile (output_file); // Open save file
  if (!savefile.is_open()) {
    cout << "Unable to open file: " << output_file << endl; // Exit if save file has not opened
    return EXIT_FAILURE;
  }
  savefile << std::setprecision(16); // Set the precision of the output to that of a double
  // Write a header for the save file
  savefile << dt << "," << T << "," << Tsave << "," << N << endl;

  // -----------------------------------------------------------  
  // Initialise variables
  double epsilon = 0.0103; // In eV
  double sigma = 3.345; // In angstroms (Å)
  double mass = 39.948; // In amu
  int save_counter = 0;
  double momentum = 0.0;
  
  update_acc(system, epsilon, sigma, mass); // Calculate acceleration at t = 0
  save_data(savefile, system); // Calculate and save data for timestep t = 0.
    
  for (int t=1; t<=T; t++) // Loop steps until it reaches total number for time steps 'T'.
  {
    std::vector<body> oldsystem = system; // Saves a copy of the current system.

    // Runs Vel_verlet integration steps to update system.
    update_vel(system, dt);
    update_pos(system, dt, boundry);
    update_acc(system, epsilon, sigma, mass);
    update_vel(system, dt); 
    
    // Pressure Calculations
    for (int i = 0; i < system.size(); i++)
    { // Determines if the particle has crossed over the plane x = 0 between timesteps.
      if (oldsystem[i].get_pos().getx() > 0 && system[i].get_pos().getx() < 0)
      {
        momentum += abs((mass * system[i].get_vel().getx())); // Summs the momentum in the x axis
      } // Makes sure the value is always positive.
    }
    
    save_counter++; // Adds 1 to Tsave step counter.

    if (save_counter == Tsave) // Checks if step counter is equal to Tsave.
    {
      save_data(savefile, system); // Saves position and velocity of each particle to the file
      save_counter = 0; // Reset Tsave counter to 0.    
    }
  }
  savefile.close();
  cout << "Average pressure of System = " << 3 * ((momentum)/(pow((boundry*2), 2) * dt * T)) << 
  " Pascals" << endl;
  cout << "End Temperature of System = " << calc_temp(system, mass) << " K" << endl;
  return EXIT_SUCCESS; 
}

// ######################################
// ###### Function implementations ######
// ######################################

void update_pos(std::vector<body> &system, double dt, double boundry)
{ 
  for (int i = 0; i < system.size(); i++)
  {
    vec pos = system[i].get_pos(); // Initialise variables for this particle.
    vec velocity = system[i].get_vel();
    
    pos += (velocity * (dt)); // Calculates new position based on new velocity.
    system[i].set_pos(pos); // Sets new position.

    // Checking if particle positions are outside the boundry. 
    if (pos.getx() >= boundry || pos.getx() <= -boundry)
    {
      system[i].set_vel(velocity.negx()); // Negates the x value of the velocity.
    }
    if (pos.gety() >= boundry || pos.gety() <= -boundry) 
    {
      system[i].set_vel(velocity.negy()); // Negates the y value of the velocity.
    }
    if (pos.getz() >= boundry || pos.getz() <= -boundry) 
    {
      system[i].set_vel(velocity.negz());// Negates the z value of the velocity.
    }
  }
}

void update_vel(std::vector<body> &system, double dt)
{ 
  for (int i = 0; i < system.size(); i++)
  {
    // Initialise variables for this planet.
    vec velocity = system[i].get_vel();
    vec acceleration = system[i].get_acc();
    
    velocity +=(acceleration * (dt / 2)); // Calculates new velocity over half a timestep.
    system[i].set_vel(velocity); // Sets new velocity to particle. In angstrom per second
  }
}

void update_acc(std::vector<body> &system, double epsilon, double sigma, double mass)
{
  for (int i = 0; i < system.size(); i++) // Loop through particles.
    {
      // Initialise variables for first particle.
      vec tot_force(0.,0.,0.);
      vec posi = system[i].get_pos();
      
      for (int j=0; j<system.size(); j++) // Loop through other particles.
      {
        if (i != j) // Particles must be different!
        {
          vec posj = system[j].get_pos(); // Initialise variables for second particle.
          
          vec difference = posi - (posj); // Calculates the vector distance between them.
          
          double length = difference.length(); // Calculates the length of the vector.
          
          // Calculates the scalar force component between the 2 particles.
          // Uses the differential equation of the Lennard-Jones potential.
          vec force = difference / (length) * (24 * (epsilon / length) * (2 * pow((sigma / length), 12) - pow((sigma / length), 6)));
          
          tot_force += (force); // Sums the different forces on the particle.
        }
      }
      // Acceleration in angstrom per unit time squared, factored down by 1.6e-19
      system[i].set_acc(tot_force / (mass)); // Calculates and sets acceleration variable of particle.
      
    }
}

double calc_temp(std::vector<body> &system, double mass)
{
  double tot_energy = 0;
  for (int i = 0; i < system.size(); i++) // Loop through particles.
  { // Factors the changes from; eV's to joules, angstroms to meters
    tot_energy += (1.6e-19 * 0.5 * mass * system[i].get_vel().lengthsq());
  }
  return(tot_energy / (1.5 * system.size() * 1.38e-23)); // Returns temperature of whole system in kelvin.
}

void read_init(std::string input_file, std::vector<body> &system)
{
  std::string line; // Declare a string to store each line
  double x, y, z, vx, vy, vz; // Doubles to store vector components
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
        read_vector3d(data_line, x, y, z); // Read the 3 components of the vector on the line
        break;
      case 2:
        read_vector3d(data_line, vx, vy, vz); // Read the 3 components of the vector on the line
        break;                              
      }
      if (line_cnt==2) // Data for one body has been extracted
      {
        line_cnt = 0; // Reset line counter
        body b(vec(x,y,z),vec(vx,vy,vz)); // Package data into body
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

void save_data(std::ofstream& savefile, const std::vector<body> &system)
{
    // Function for saving the simulation data to file.
    // Loop over the bodies:
    for(int p = 0; p < system.size(); p++)
    { 
      // Output position and velocity for each body:
      savefile << system[p].get_pos() << "," << system[p].get_vel() << endl;
    }
}
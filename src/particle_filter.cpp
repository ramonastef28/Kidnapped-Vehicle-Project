/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include  <math.h>
#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
        num_particles = 400;
        // noise generation
        default_random_engine gen;
        normal_distribution<double> x_part(x, std[0]);
        normal_distribution<double> y_part(y, std[1]);
        normal_distribution<double> theta_part(theta, std[2]);
  
	for (unsigned int i=0; i< num_particles; ++i)
		{
		Particle particle;
		particle.id = i;
		particle.x = x_part(gen);
		particle.y = y_part(gen);
		particle.theta = theta_part(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
		weights.push_back(particle.weight);
		//cout << "Particle " << particle.id << "x " << particle.x << "y " << particle.y << "theta " << particle.theta << endl;
		}
	is_initialized = true;
}	

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	for (unsigned int i=0; i<num_particles; ++i)
		{
		double theta_new= particles[i].theta + yaw_rate*delta_t;
		if (abs(yaw_rate) > 0.00001) {
			particles[i].x = particles[i].x + (velocity/yaw_rate) * (sin(theta_new) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(theta_new));
			}
		else{
			particles[i].x = velocity*delta_t*cos(particles[i].theta);
			particles[i].y = velocity*delta_t*sin(particles[i].theta);
			}	
		particles[i].theta = theta_new;
		//cout << "particle before" << particles[i].x << endl;
		//adding Gaussian random noise to each particle
		normal_distribution<double> particle_x(particles[i].x, std_pos[0]);
		normal_distribution<double> particle_y(particles[i].y, std_pos[1]);
		normal_distribution<double> particle_theta(particles[i].theta , std_pos[2]);

		particles[i].x =  particle_x(gen);
		particles[i].y =  particle_y(gen);
		particles[i].theta = particle_theta(gen);	
		}
}

//void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, Map map_landmarks, std::vector<LandmarkObs>& LandmarkAssociation){ 
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
     
      	//cout << "size global obs" << predicted.size() << endl;
      	//cout << " size map_landmarks" << map_landmarks.landmark_list.size() << endl;  	
	//for (int i=0; i < predicted.size(); i++){
	//	cout << "x global" << predicted[i].x << endl;}
	//for (unsigned int i=0; i < map_landmarks.landmark_list.size(); i++){
        //	 cout << "landmark x " << i << " " << map_landmarks.landmark_list[i].x_f << endl;}		


	for (unsigned int i=0; i < predicted.size(); ++i){
		LandmarkObs LandmarkA;
		double min_dist;

		LandmarkA.id  = map_landmarks.landmark_list[0].id_i;
		LandmarkA.x = map_landmarks.landmark_list[0].x_f;
		LandmarkA.y = map_landmarks.landmark_list[0].y_f;
		min_dist = dist(LandmarkA.x, LandmarkA.y, predicted[i].x, predicted[i].y);
		for (unsigned int j=1; j <  map_landmarks.landmark_list.size(); ++j) {
			double dist_land = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f,
			predicted[i].x, predicted[i].y);
			if (dist_land < min_dist) {
				min_dist = dist_land;
				LandmarkA.id = map_landmarks.landmark_list[j].id_i;
				LandmarkA.x = map_landmarks.landmark_list[j].x_f;
				LandmarkA.y = map_landmarks.landmark_list[j].y_f;
			}
		}
		LandmarkAssociation.push_back(LandmarkA);
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	
	//cout << "size observations " << observations.size() << endl;
	//cout << "size map " << map_landmarks.landmark_list.size() << endl;

	//for (unsigned int i=0; i < observations.size(); i++){
	//	cout << "obs x " << observations[i].x << endl; }  

	//for (unsigned int i=0; i < map_landmarks.landmark_list.size(); i++){
	//	cout << "landmark x " << i << "" << map_landmarks.landmark_list[i].x_f << endl;}

        for (unsigned int i=0; i<num_particles; i++){
		std::vector<LandmarkObs> globalObs;
		std::vector<LandmarkObs> LandmarkAssociation;
        	for (unsigned int j=0; j< observations.size(); ++j){
			LandmarkObs tempObs;
			tempObs.id = observations[j].id;
			tempObs.x = observations[j].x*cos(particles[i].theta) - 
				observations[j].y *sin(particles[i].theta) + particles[i].x;
			tempObs.y = observations[j].x*sin(particles[i].theta) + 
				observations[j].y*cos(particles[i].theta) + particles[i].y;
			globalObs.push_back(tempObs);	
		}
	     	dataAssociation(globalObs, map_landmarks, LandmarkAssociation);	        
				
		double weight_final = 1;
		double sigma_x = std_landmark[0];
		double sigma_y = std_landmark[1];
		for (int j =0; j < globalObs.size(); ++j){
			double delta_x = globalObs[j].x - LandmarkAssociation[j].x;
			double delta_y = globalObs[j].y - LandmarkAssociation[j].y;
			double prob1 = 1.0/(2 * M_PI * sigma_x * sigma_y);
			double prob2 = pow(delta_x, 2)/pow(sigma_x, 2) + pow(delta_y,2)/pow(sigma_y, 2);
			prob2 = exp(-0.5*prob2);
			weight_final = weight_final*prob1*prob2;
		}
		particles[i].weight = weight_final;
		weights[i] = weight_final;		
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<> dd(weights.begin(), weights.end());
	vector<Particle> newParticles;
	for (int i = 0; i < num_particles; ++i){
		newParticles.push_back(particles[dd(gen)]);
	}	
	particles = newParticles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

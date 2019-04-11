/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */


#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  // Random number generator
  std::default_random_engine gen;
  
  // Normal distributions for initial coordinates
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  // Initialize the weights and particles
  particles.resize(num_particles);
  weights.resize(num_particles);
  for (int i = 0; i < num_particles; i++)
  {
    double x = dist_x(gen);
    double y = dist_y(gen);
    double theta = dist_theta(gen);
    double weight = 1.;
    
    Particle particle;
    particle.id = i;
    particle.x = x;
    particle.y = y;
    particle.theta = theta;    
    particle.weight = weight;
    
    particles[i] = particle;
    weights[i] = weight;
  }
  
  is_initialized = true;  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0.0, std_pos[0]);
  std::normal_distribution<double> dist_y(0.0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0.0, std_pos[2]);  

  	for (int i=0; i < num_particles; i++) {
      
      if (fabs(yaw_rate) < 0.000001) {
			particles[i].x += cos(particles[i].theta) * velocity * delta_t;
			particles[i].y += sin(particles[i].theta) * velocity * delta_t;
			// particles[i].theta unchanged if yaw_rate is too small
		} else {
			double delta_theta = yaw_rate * delta_t;
			particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + delta_theta) - sin(particles[i].theta));
			particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + delta_theta));
			particles[i].theta += yaw_rate * delta_t;
		}
		// Add noise
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].y += dist_theta(gen);      
    }  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

	for (int i = 0; i < observations.size(); i++) {
      	LandmarkObs& observation = observations[i];
      
		double dist_to_predict = -1; 
		int predict_id = -1;
		for (int j = 0; j < predicted.size(); j++) {
          	LandmarkObs predict = predicted[j];
          
			double distance = dist(observation.x, observation.y, predict.x, predict.y);
			if (dist_to_predict == -1 || distance < dist_to_predict) {
				dist_to_predict = distance;
				predict_id = predict.id;
			}
		}
		observation.id = predict_id;
	}    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

	for (int i = 0; i < num_particles; i++) {
      
		double particle_x = particles[i].x;
		double particle_y = particles[i].y;
		double particle_theta = particles[i].theta;
      
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		particles[i].associations.clear();
      
		vector<LandmarkObs> predicted;
		for (Map::single_landmark_s landmark : map_landmarks.landmark_list) {
			// assuming sensor range squarely in both x and y direction
			if (fabs(landmark.x_f-particle_x) <= sensor_range && fabs(landmark.y_f-particle_y) <= sensor_range) {
				LandmarkObs predict;
				predict.id = landmark.id_i;
				predict.x = landmark.x_f;
				predict.y = landmark.y_f;
				predicted.push_back(predict);
			}
		}
      
		// list of actual observations, transformed into map coordinates
		vector<LandmarkObs> observations_transformed;
		for (LandmarkObs observation: observations) {
			LandmarkObs observation_transformed;
			observation_transformed.x = cos(particle_theta)*observation.x - sin(particle_theta)*observation.y + particle_x;;
			observation_transformed.y = sin(particle_theta)*observation.x + cos(particle_theta)*observation.y + particle_y;
			observation_transformed.id = observation.id;
			observations_transformed.push_back(observation_transformed);
		}
      
		// associate predicted landmarks with observed (measurement) landmarks
		dataAssociation(predicted, observations_transformed);

		particles[i].weight = 1.0;
		// find actual x and y coordinates for each observed landmark
		for (LandmarkObs observation: observations_transformed) {
          
			double predict_x, predict_y;
			for (LandmarkObs predict: predicted) {
				if (predict.id == observation.id) {
					predict_x = predict.x;
					predict_y = predict.y;
					break;
				}
			}
          
			// calculate new weight with multi variate Gaussian probability density function
			double sigma_x = std_landmark[0];
			double sigma_y = std_landmark[1];
			double normalizer = (1 / (2 * M_PI * sigma_x * sigma_y));
			double exponent = pow(predict_x-observation.x, 2) / (2 * pow(sigma_x, 2)) + 
								pow(predict_y-observation.y, 2) / (2 * pow(sigma_y, 2));
			double weight = normalizer * exp(-exponent);
          
			particles[i].weight *= weight;
			weights[i] = particles[i].weight;

			particles[i].associations.push_back(observation.id);
			particles[i].sense_x.push_back(observation.x);
			particles[i].sense_y.push_back(observation.y);
		}
	}
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    // Create a discrete distribution and a generator
    std::default_random_engine gen;
    std::discrete_distribution<int> distribution {weights.begin(), weights.end()};

    // Create vector for new particles
    std::vector<Particle> new_particles;
     // Draw new particles based on distribution
    for (int i = 0; i < num_particles; i++) {
        // Generate particle index using the distribution
        int cur_particle_i = distribution(gen);

        // Get a particle from the list of particles
        Particle cur_particle = particles[cur_particle_i];

        // Push back particle to the list of new particles
        new_particles.push_back(cur_particle);
    }

    // Make new particles current particles
    particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
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
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  cout << "\n--------------- Initializing particle filter ---------------" << endl;

  // Particle filter size
  num_particles = 250;

  // GPS/measurement uncertainties
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // Gaussian generator
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  for (int i = 0; i < num_particles; i++) {
    // Generate gaussian samples
    double sample_x = dist_x(gen);
    double sample_y = dist_y(gen);
    double sample_theta = dist_theta(gen);

    Particle particle;
    particle.id = i;
    particle.x = sample_x;
    particle.y = sample_y;
    particle.theta = sample_theta;
    particle.weight = 1.0f;

    weights.push_back(particle.weight);
    particles.push_back(particle);
  }
  cout << num_particles << " particles created" << endl;

  is_initialized = true;

  cout << "Particle filter initialization done" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // Avoid division by 0
  if (fabs(yaw_rate < 0.00001)) {
    if (yaw_rate < 0) {
      yaw_rate = -0.00001;
    } else {
      yaw_rate = 0.00001;
    }
  }

  // GPS/measurement uncertainties
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  // Generate gaussian samples
  normal_distribution<double> noise_x(0, std_x);
  normal_distribution<double> noise_y(0, std_y);
  normal_distribution<double> noise_theta(0, std_theta);

  for (int i = 0; i < num_particles; i++) {
    Particle &particle = particles[i];
    // Make prediction
    double x_t0 = particle.x;
    double y_t0 = particle.y;
    double theta_t0 = particle.theta;

    // cout << "before particle update:" << endl;
    // cout << "x, y, theta=" << x_t0 << "," << y_t0 << "," << theta_t0 << endl;
    // cout << "noise x, y, theta=" << noise_x(gen) << "," << noise_y(gen) << "," << noise_y(gen) << endl;

    double x_t1_with_noise =
        x_t0 + velocity / yaw_rate * (sin(theta_t0 + yaw_rate * delta_t) - sin(theta_t0)) + noise_x(gen);
    double y_t1_with_noise =
        y_t0 + velocity / yaw_rate * (-cos(theta_t0 + yaw_rate * delta_t) + cos(theta_t0)) + noise_y(gen);
    double theta_t1_with_noise = theta_t0 + yaw_rate * delta_t + noise_theta(gen);

    // Update particle filter
    particle.x = x_t1_with_noise;
    particle.y = y_t1_with_noise;
    particle.theta = theta_t1_with_noise;

    // cout << "after particle update:" << endl;
    // cout << "x, y, theta=" << particle.x << "," << particle.y << "," << particle.theta << endl;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html
  double std_landmark_x = std_landmark[0];
  double std_landmark_y = std_landmark[1];

  double var_landmark_x = std_landmark_x * std_landmark_x;
  double var_landmark_y = std_landmark_y * std_landmark_y;

  double covar_landmark_xy = std_landmark_x * std_landmark_y;

  // Translate observation to map coord system
  long double total_weight = 0;
  for (int i = 0; i < num_particles; i++) {
    Particle &particle = particles[i];

    long double weight = 1.0;
    // cout << observations.size() << " observations..." << endl;

    for (int j = 0; j < observations.size(); j++) {
      LandmarkObs obs = observations[j];

      // Get predictions of landmarks in map coord system
      double pred_x = particle.x + cos(particle.theta) * obs.x - sin(particle.theta) * obs.y;
      double pred_y = particle.y + sin(particle.theta) * obs.x + cos(particle.theta) * obs.y;

      Map::single_landmark_s nearest_landmark;
      double min_dist = sensor_range;
      double d = 0;

      // Find out which is the nearest landmark
      for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
        Map::single_landmark_s landmark = map_landmarks.landmark_list[k];

        // Get the distance to prediction
        d = dist(pred_x, pred_y, landmark.x_f, landmark.y_f);
        if (d < min_dist) {
          min_dist = d;
          nearest_landmark = landmark;
        }
      }

      // cout << "nearest landmark id: " << nearest_landmark.id_i << endl;

      // Get differnce between prediction and matched landmark
      double x_diff = pred_x - nearest_landmark.x_f;
      double y_diff = pred_y - nearest_landmark.y_f;

      // Observation weight
      weight = (1 / (2 * M_PI * covar_landmark_xy)) *
               exp(-(pow(x_diff, 2) / (2 * var_landmark_x) + (pow(y_diff, 2) / (2 * var_landmark_y))));
//      cout << "x_diff:" << x_diff << " y_diff:" << y_diff << " std_x:" << std_landmark_x << " std_y:" << std_landmark_y
//           << " weight:" << weight << endl;

      // Total prob
      weight *= weight;
    }
    // Assign weight to particle
    // cout << "particle final weight:" << weight << endl;
    particle.weight = double(weight);

    // Add up to the total weight
    total_weight += weight;

  }

  // Normalize particle weights
//  cout << "total weight of all particles:" << total_weight << endl;
//  for (int i = 0; i < num_particles; i++) {
//    // Update particle weight
//    particles[i].weight /= total_weight;
//    // Update weight in PF array
//    weights[i] /= total_weight;
//  }

}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Use a discrete distribution where probability of each random idx being selected is propotional to its weights
  discrete_distribution<> dd_by_weights(weights.begin(), weights.end());

  // Create new particle list
  vector<Particle> new_particles;
  // Sample particles from current particle list and add to new list
  for (int i = 0; i < num_particles; i++)
    new_particles.push_back(particles[dd_by_weights(gen)]);

  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
                                         const std::vector<double> &sense_x, const std::vector<double> &sense_y) {
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

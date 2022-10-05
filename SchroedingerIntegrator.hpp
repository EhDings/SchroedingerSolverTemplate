//
// Created by Michael Eichberg on 05.10.22.
//

#include <functional>
#include <utility>
#include <vector>

#include <boost/numeric/odeint.hpp>

#ifndef _SCHROEDINGERINTEGRATOR_HPP_
#define _SCHROEDINGERINTEGRATOR_HPP_

using boost::numeric::odeint::integrate_const;
using boost::numeric::odeint::make_controlled;
using boost::numeric::odeint::runge_kutta_fehlberg78;

class SchroedingerIntegrator {
 public:
  SchroedingerIntegrator(double energy, double reduced_mass, int L, std::function<double(double)> potential, double rmin, size_t tracked_index)
  : m_energy{energy}, m_reduced_mass{reduced_mass}, m_L{L}, m_potential{std::move(potential)}, m_last_position{rmin}, m_tracked_index{tracked_index}
  {
    /*
     * Initialize here m_last_state at rmin.
     */
  }

  void do_integrate_range(double range, double stepsize, double abs_error, double rel_error, bool do_track)
  {
    /*
     * range: Integrate from m_last_position to m_last_position+range.
     * stepsize: Initial stepsize and distance between tracked positions.
     * abs, rel_error: Tolerances for adaptive stepsize.
     * do_track: If true, m_positions and m_trajectory will save all points.
     *           If false, only information is saved in m_last_position and m_last_state.
     */
    if (do_track)
      integrate_const (make_controlled<runge_kutta_fehlberg78<std::vector<double>>> (abs_error, rel_error), m_state_velocity, m_last_state, m_last_position,
        m_last_position + range, stepsize, m_tracker);
    else
    integrate_const (make_controlled<runge_kutta_fehlberg78<std::vector<double>>> (abs_error, rel_error), m_state_velocity, m_last_state, m_last_position,
        m_last_position + range, range, m_tracker);
    m_last_position = m_positions.back();
  }

  double get_last_point() { return m_last_state[m_tracked_index]; }

 private:
  class TrajectoryTracker{
   public:
    std::vector<double> &m_positions;
    std::vector<double> &m_trajectory;
    size_t &m_tracked_index;

    TrajectoryTracker (std::vector<double> &positions, std::vector<double> &trajectory, size_t &tracked_index)
        : m_positions{positions}, m_trajectory{trajectory}, m_tracked_index{tracked_index}
    {}

    void operator() (const std::vector<double> &state, double pos) const
    {
      /*
       * Include routine to save relevant information in positions and trajectory vectors here.
       */
    }
  };
  class StateVelocity{
   public:
    double &m_energy;
    double &m_reduced_mass;
    int &m_L;
    std::function<double (double)> &m_potential;

    explicit StateVelocity (double &energy, double &reduced_mass, int &L, std::function<double (double)> &potential)
        : m_energy{energy}, m_reduced_mass{reduced_mass}, m_L{L}, m_potential{potential}
    {}

    void operator() (const std::vector<double> &state, std::vector<double> &ODE_vec, double pos) const
    {
      /*
       * state can be e.g. (psi, psi', psi'').
       * ODE_vec are the ODEs for psi, psi', psi'', which will be defined here.
       */
    }
  };
  double m_energy{};
  double m_reduced_mass{}; // Reduced mass in the kinetic operator.
  int m_L{}; // Angular momentum.
  std::function<double(double)> m_potential{}; // Functor with operator() giving the value of V(r) at r.
  double m_last_position{}; // Last r.
  size_t m_tracked_index{}; // The index of the quantity relevant for boundary conditions.
  std::vector<double> m_last_state{}; // Last state of wave function vector (psi, psi', psi'').
  std::vector<double> m_positions{}; // r's at which wave function has been evaluated.
  std::vector<double> m_trajectory{}; // Values of the quantity relevant for boundary conditions at all r's.
  StateVelocity m_state_velocity{m_energy, m_reduced_mass, m_L, m_potential};
  TrajectoryTracker m_tracker{m_positions, m_trajectory, m_tracked_index};
};

#endif //_SCHROEDINGERINTEGRATOR_HPP_

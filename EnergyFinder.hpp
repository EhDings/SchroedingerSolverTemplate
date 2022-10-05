//
// Created by Michael Eichberg on 05.10.22.
//

#include <functional>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>

#include <boost/math/tools/roots.hpp>

#include "SchroedingerIntegrator.hpp"

#ifndef _ENERGYFINDER_HPP_
#define _ENERGYFINDER_HPP_

using boost::math::tools::toms748_solve;
using boost::math::tools::eps_tolerance;

class EnergyFinder {
 public:
  EnergyFinder ()
  {}

  double do_scan (double energy)
  {
    return m_root_criterium (energy);
  }

  void do_collect_brackets (double start, double stepsize, size_t max_n, bool check_bottom, size_t max_scans = 1000)
  {
    if (stepsize == 0.0)
      throw std::runtime_error ("Error: Scan with stepsize = 0");
    size_t i1{};
    size_t i2{};
    double old_energy;
    double energy{start};
    double old_value;
    double new_value{do_scan (energy)};
    if (check_bottom)
      {
        double gap_size{};
        double last_gap_size{fabs (stepsize) * (double) max_scans};
        while (true)
          {
            gap_size += fabs (stepsize);
            if (gap_size > last_gap_size * 2.0)
              {
                std::cout << "We think we've reached the bottom." << "\n";
                break;
              }
            old_energy = energy;
            energy = start - fabs (stepsize) * (double) (i1 + 1);
            old_value = new_value;
            new_value = do_scan (energy);
            if (new_value * old_value < 0.0)
              {
                i2++;
                m_brackets[energy] = old_energy;
                if (m_brackets.size () > 2)
                  {
                    std::map<double, double>::iterator bracket_iterator{m_brackets.begin ()};
                    bracket_iterator++;
                    last_gap_size = fabs (bracket_iterator->second - m_brackets.begin ()->first);
                  }
              }
          }
      }
    energy = start;
    new_value = do_scan (energy);
    while (i2 < max_n)
      {
        old_energy = energy;
        energy = start + stepsize * (double) (i1 + 1);
        old_value = new_value;
        new_value = do_scan (energy);
        i1++;
        if (new_value * old_value < 0.0)
          {
            i2++;
            if (stepsize > 0)
              m_brackets[old_energy] = energy;
            else
              m_brackets[energy] = old_energy;
          }
        if (i1 > max_scans)
          {
            std::cerr << "Warning: Maximum number of scans reached, #max_scans = " << max_scans << "\n";
            break;
          }
      }
  }

  double do_shoot (double minimum, double maximum)
  {
    if (minimum > maximum)
      {
        double tmp{minimum};
        minimum = maximum;
        maximum = tmp;
      }
    const std::uintmax_t max_it{20};
    std::uintmax_t it{max_it};
    int digits{std::numeric_limits<double>::digits};
    int get_digits{digits - 3};
    eps_tolerance<double> tol (get_digits);
    std::cout << "Find root in bracket ( minimum, maximum ) = ( " << minimum << ",\t" << maximum << " )\n";
    std::pair<double, double> bracket{toms748_solve (m_root_criterium, minimum, maximum, tol, it)};
    return (bracket.first + bracket.second) / 2.0;
  }

  void do_shoot_brackets ()
  {
    for (const auto &bracket: m_brackets)
      {
        m_energies.push_back (do_shoot (bracket.first, bracket.second));
      }
  }

 private:
  class RootCriterium {
    /*
     * Nested class to define the quantity of which we want to find the root.
     */
   public:
    int &m_L;
    size_t &m_tracked_index;
    double &m_rmin;
    double &m_reduced_mass;
    double &m_max_range;
    double &m_stepsize;
    double &m_abs_error;
    double &m_rel_error;
    std::function<double (double)> &m_potential;

    RootCriterium (std::function<double (double)> &potential, double &rmin, double &reduced_mass, int &L, size_t &tracked_index, double &max_range, double &stepsize, double &abs_error, double &rel_error)
        : m_L{L}, m_tracked_index{tracked_index}, m_rmin{rmin}, m_reduced_mass{reduced_mass}, m_max_range{
        max_range}, m_stepsize{stepsize}, m_abs_error{abs_error}, m_rel_error{rel_error}, m_potential{potential}
    {}

    double operator() (double energy) const
    {
      /*
       * Here, we want to find the root of the wave function psi at rmax with respect to the energy.
       */
      size_t i1{};
      while (energy > m_potential (m_max_range + m_rmin))
        {
          i1++;
          m_max_range += 1.0;
          if (i1 > 100)
            break;
        }
      if (m_max_range > 5.0)
        {
          if (energy > m_potential (m_max_range + m_rmin - 5.0))
            m_max_range += 10.0;
        }
      SchroedingerIntegrator schroedinger_integrator (energy, m_reduced_mass, m_L, m_potential, m_rmin, m_tracked_index);
      schroedinger_integrator.do_integrate_range (m_max_range, m_stepsize, m_abs_error, m_rel_error, false);
      return schroedinger_integrator.get_last_point ();
    }
  };
  int m_L{};
  size_t m_tracked_index{};
  double m_rmin{};
  double m_reduced_mass{};
  double m_max_range{};
  double m_stepsize{};
  double m_abs_error{};
  double m_rel_error{};
  std::map<double, double> m_brackets{};
  std::vector<double> m_energies{};
  std::function<double (double)> m_potential{};
  RootCriterium m_root_criterium{m_potential, m_rmin, m_reduced_mass, m_L, m_tracked_index, m_max_range, m_stepsize,
                                 m_abs_error, m_rel_error};
};

#endif //_ENERGYFINDER_HPP_

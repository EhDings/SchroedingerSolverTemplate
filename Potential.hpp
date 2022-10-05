//
// Created by Michael Eichberg on 05.10.22.
//

/*
 * Simple class for defining a potential functor with specific parameters.
 */

#ifndef _POTENTIAL_HPP_
#define _POTENTIAL_HPP_

class Potential {
 public:
  Potential(double e, double sigma) : m_e{e}, m_sigma{sigma} {} // Initialize the parameters of the potential.

  double operator() (double pos) { return - m_e / pos + m_sigma * pos; } // Makes this class a functor.

  Potential &set_param1(double e)
  {
    m_e = e;
    return *this;
  }

  Potential &set_param2(double sigma)
  {
    m_sigma = sigma;
    return *this;
  }

  double get_e() { return m_e; }
  double get_sigma() { return m_sigma; }

 private:
  double m_e{};
  double m_sigma{};
};

#endif //_POTENTIAL_HPP_

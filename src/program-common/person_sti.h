#ifndef PERSON_STI_H
#define PERSON_STI_H

#include <assert.h>

class Person;
class ConfigSettings;
class ConfigWriter;
class GslRandomNumberGenerator;
class ProbabilityDistribution;


class Person_STI {
public:
  enum InfectionType { None, Partner, Seed };
  enum InfectionSite { Vaginal, Rectal, Urethral};
  
  Person_STI(Person *pSelf): m_pSelf(pSelf), m_infectionTime(-1e200), m_pInfectionOrigin(0), m_infectionType(None), m_infectionSite(Vaginal) {}
  
  virtual ~Person_STI() {}
  
  virtual bool isInfected() const = 0;
  virtual bool isInfectious() const = 0;
  
  InfectionType getInfectionType() const { return m_infectionType; }
  InfectionSite getInfectionSite() const { return m_infectionSite; }
  double getInfectionTime() const { assert(isInfected()); return m_infectionTime; }
  Person *getInfectionOrigin() const { assert(isInfected()); return m_pInfectionOrigin; }
  double getRecoveryTime() const                          { return m_recoveryTime; }
  bool isTreated() const                                  { return m_treated; }
  
  virtual void setInfected(double t, Person *pOrigin, InfectionType iType, InfectionSite iSite) = 0;
  virtual void progress(double t, bool treatInd) = 0;
  
  
  
protected:
  const Person *m_pSelf;
  
  double m_infectionTime;
  double m_recoveryTime;
  bool m_treated;
  Person *m_pInfectionOrigin;
  InfectionType m_infectionType;
  InfectionSite m_infectionSite;
};

#endif // PERSON_STI_H
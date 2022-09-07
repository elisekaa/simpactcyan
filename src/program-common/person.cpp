#include "person.h"
#include "personimpl.h"
#include "configsettings.h"
#include "configwriter.h"
#include "configdistributionhelper.h"
#include "debugwarning.h"
#include "logsystem.h"
#include "simpactevent.h"
#include "discretedistribution2d.h"
#include "jsonconfig.h"
#include "configfunctions.h"
#include <stdlib.h>
#include <limits>

using namespace std;

Person::Person(double dateOfBirth, Gender g) : PersonBase(g, dateOfBirth), m_relations(this), m_hiv(this),
	                                           m_hsv2(this), m_gonorrhea(this), m_condom_use_probability_concordant(0),
											   m_condom_use_probability_discordant(0), m_health_seeking_propensity(0)
{
	assert(g == Male || g == Female);

	assert(m_pPopDist);
	assert(m_pHealthSeekingPropensityDist);
	assert(m_pCondomUseProbDist);

	Point2D loc = m_pPopDist->pickPoint();
	assert(loc.x == loc.x && loc.y == loc.y); // check for NaN
	setLocation(loc, 0);

	double hsp = m_pHealthSeekingPropensityDist->pickNumber();
	assert(hsp == hsp); // check for NaN
	setHealthSeekingPropensity(hsp);

	Point2D condomUseProbs = m_pCondomUseProbDist->pickPoint();
	assert(condomUseProbs.x == condomUseProbs.x && condomUseProbs.y == condomUseProbs.y); // check for NaN
	m_condom_use_probability_concordant = condomUseProbs.x;
	m_condom_use_probability_discordant = condomUseProbs.y; // TODO use setters?

	m_pPersonImpl = new PersonImpl(*this);
}

Person::~Person()
{
	delete m_pPersonImpl;
}

ProbabilityDistribution2D *Person::m_pPopDist = 0;
double Person::m_popDistWidth = 0;
double Person::m_popDistHeight = 0;

ProbabilityDistribution *Person::m_pHealthSeekingPropensityDist = 0;

ProbabilityDistribution2D *Person::m_pCondomUseProbDist = 0;

double Person::getCondomUseProbability(bool isPartnerDiagnosed) const
{
	bool amIDiagnosed = m_hiv.isDiagnosed();

	if (isPartnerDiagnosed == amIDiagnosed) {
		// Concordant HIV status
		return m_condom_use_probability_concordant;
	} else  {
		// Discordant HIV status
		return m_condom_use_probability_discordant;
	}
	// TODO base on ART use?
	// TODO base on PreP use?
	// TODO concordant + and concordant - different?
	// FIXME should this be framed as 'probability' in the context of hazard functions?
}

bool Person::isInfectedWithSTI() const
{
	return (m_gonorrhea.isInfected(), m_hsv2.isInfected());
}

void Person::processConfig(ConfigSettings &config, GslRandomNumberGenerator *pRndGen)
{
	assert(pRndGen != 0);

	// Population distribution
	delete m_pPopDist;
	m_pPopDist = getDistribution2DFromConfig(config, pRndGen, "person.geo");

	// Health-seeking propensity distribution
	delete m_pHealthSeekingPropensityDist;
	m_pHealthSeekingPropensityDist = getDistributionFromConfig(config, pRndGen, "person.healthseekingpropensity");

	// Condom use probability distributions
	delete m_pCondomUseProbDist;
	m_pCondomUseProbDist = getDistribution2DFromConfig(config, pRndGen, "person.condomuse");
}

void Person::obtainConfig(ConfigWriter &config)
{
	assert(m_pPopDist);
	addDistribution2DToConfig(m_pPopDist, config, "person.geo");
	assert(m_pHealthSeekingPropensityDist);
	addDistributionToConfig(m_pHealthSeekingPropensityDist, config, "person.healthseekingpropensity");
	assert(m_pCondomUseProbDist);
	addDistribution2DToConfig(m_pCondomUseProbDist, config, "person.condomuse");
}

void Person::writeToPersonLog()
{
	double infinity = numeric_limits<double>::infinity();
	double NaN = numeric_limits<double>::quiet_NaN();

	int id = (int)getPersonID(); // TODO: should fit in an 'int' (easier for output)
	int gender = (isMan())?0:1;
	double timeOfBirth = getDateOfBirth();
	double timeOfDeath = (hasDied())?getTimeOfDeath():infinity;

	Man *pFather = getFather();
	Woman *pMother = getMother();
	int fatherID = (pFather != 0) ? (int)pFather->getPersonID() : (-1); // TODO: cast should be ok
	int motherID = (pMother != 0) ? (int)pMother->getPersonID() : (-1);

	// TODO: Currently not keeping track of children

	double debutTime = (isSexuallyActive())? m_relations.getDebutTime():infinity;
	double formationEagerness = getFormationEagernessParameter();
	double formationEagernessMSM = (isMan())?getFormationEagernessParameterMSM():NaN;
	
	double infectionTime = (m_hiv.isInfected())? m_hiv.getInfectionTime() : infinity;
	Person *pOrigin = (m_hiv.isInfected())? m_hiv.getInfectionOrigin() : 0;
	int origin = (pOrigin != 0) ? (int)pOrigin->getPersonID() : (-1); // TODO: cast should be ok
	
	int infectionType = 0;
	switch(m_hiv.getInfectionType())
	{
	case Person_HIV::None:
		infectionType = -1;
		break;
	case Person_HIV::Seed:
		infectionType = 0;
		break;
	case Person_HIV::Partner:
		infectionType = 1;
		break;
	case Person_HIV::Mother:
		infectionType = 2;
		break;
	default: // Unknown, but don't abort the program at this point
		infectionType = 10000 + (int)m_hiv.getInfectionType();
	}

	double log10SPVLoriginal = (m_hiv.isInfected()) ? std::log10(m_hiv.getOriginalViralLoad()) : -infinity;
	double treatmentTime = (m_hiv.isInfected() && m_hiv.hasLoweredViralLoad()) ? m_hiv.getLastTreatmentStartTime() : infinity;

	int aidsDeath = -1;
	if (hasDied())
	{
		aidsDeath = 0;
		if (m_hiv.wasAIDSDeath())
			aidsDeath = 1;
	}

	double hsv2InfectionTime = (m_hsv2.isInfected())? m_hsv2.getInfectionTime() : infinity;
	Person *pHSV2Origin = (m_hsv2.isInfected()) ? m_hsv2.getInfectionOrigin() : 0;
	int hsv2origin = (pHSV2Origin != 0) ? (int)pHSV2Origin->getPersonID() : (-1); // TODO: cast should be ok

	double cd4AtInfection = (m_hiv.isInfected())?m_hiv.getCD4CountAtInfectionStart() : (-1);
	double cd4AtDeath = (m_hiv.isInfected())?m_hiv.getCD4CountAtDeath() : (-1);

	LogPerson.print("%d,%d,%10.10f,%10.10f,%d,%d,%10.10f,%10.10f,%10.10f,%10.10f,%d,%d,%10.10f,%10.10f,%10.10f,%10.10f,%d,%10.10f,%d,%10.10f,%10.10f,%10.10f",
		        id, gender, timeOfBirth, timeOfDeath, fatherID, motherID, debutTime,
		        formationEagerness,formationEagernessMSM,
		        infectionTime, origin, infectionType, log10SPVLoriginal, treatmentTime,
				m_location.x, m_location.y, aidsDeath,
				hsv2InfectionTime, hsv2origin, cd4AtInfection, cd4AtDeath, m_health_seeking_propensity);
}

void Person::writeToLocationLog(double tNow)
{
	LogLocation.print("%10.10f,%d,%10.10f,%10.10f", tNow, (int)getPersonID(), m_location.x, m_location.y);
}

void Person::writeToTreatmentLog(double dropoutTime, bool justDied)
{
	int id = (int)getPersonID(); // TODO: should fit in an 'int' (easier for output)
	int gender = (isMan())?0:1;
	int justDiedInt = (justDied)?1:0;
	double lastTreatmentStartTime = m_hiv.getLastTreatmentStartTime();
	double lastCD4AtDiagnosis = m_hiv.getLastCD4CountAtDiagnosis();
	double lastCD4AtTreatmentStart = m_hiv.getLastCD4CountAtARTStart();

	assert(m_hiv.hasLoweredViralLoad());
	assert(lastTreatmentStartTime >= 0);

	LogTreatment.print("%d,%d,%10.10f,%10.10f,%d,%10.10f", id, gender, lastTreatmentStartTime, 
	                                                       dropoutTime, justDiedInt, lastCD4AtDiagnosis, lastCD4AtTreatmentStart);
}

Man::Man(double dateOfBirth) : Person(dateOfBirth, Male)
{
}

Man::~Man()
{
}

Woman::Woman(double dateOfBirth) : Person(dateOfBirth, Female)
{
	m_pregnant = false;
}

Woman::~Woman()
{
}

ConfigFunctions personConfigFunctions(Person::processConfig, Person::obtainConfig, "Person");


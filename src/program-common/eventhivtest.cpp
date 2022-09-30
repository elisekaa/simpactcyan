#include "configsettings.h"
#include "configwriter.h"
#include "eventmonitoring.h"
#include "configdistributionhelper.h"
#include "gslrandomnumbergenerator.h"
#include "jsonconfig.h"
#include "configfunctions.h"
#include "util.h"
#include <iostream>

#include "eventhivtest.h"

using namespace std;

EventHIVTest::EventHIVTest(Person *pPerson, bool scheduleImmediately) : SimpactEvent(pPerson), m_scheduleImmediately(scheduleImmediately)
{
}

EventHIVTest::~EventHIVTest()
{
}

string EventHIVTest::getDescription(double tNow) const
{
	return strprintf("HIV test event for %s", getPerson(0)->getName().c_str());
}

void EventHIVTest::writeLogs(const SimpactPopulation &pop, double tNow) const
{
	Person *pPerson = getPerson(0);
	writeEventLogStart(true, "hivtest", tNow, pPerson, 0);
}

void EventHIVTest::markOtherAffectedPeople(const PopulationStateInterface &population)
{
	Person *pPerson = getPerson(0);

	// Infected partners (who possibly have a diagnosis event, of which
	// the hazard depends on the number of diagnosed partners), are also
	// affected!
	int numRel = pPerson->getNumberOfRelationships();

	pPerson->startRelationshipIteration();
	for (int i = 0 ; i < numRel ; i++)
	{
		double tDummy;
		Person *pPartner = pPerson->getNextRelationshipPartner(tDummy);

		if (pPartner->hiv().isInfected())
			population.markAffectedPerson(pPartner);
	}

#ifndef NDEBUG
	// Double check that the iteration is done
	double tDummy;
	assert(pPerson->getNextRelationshipPartner(tDummy) == 0);
#endif // NDEBUG
}

void EventHIVTest::fire(Algorithm *pAlgorithm, State *pState, double t)
{
	SimpactPopulation &population = SIMPACTPOPULATION(pState);
	Person *pPerson = getPerson(0);

	pPerson->hiv().setLastTestTime(t);

	// If person is HIV infected, diagnose & start monitoring,
	// else if person is not infected, schedule new test if not on PreP.
	// If on PreP, individuals are tested on each PreP screening.
	if (pPerson->hiv().isInfected()) {
		assert(!pPerson->hiv().hasLoweredViralLoad());

		pPerson->hiv().increaseDiagnoseCount(t);

		// If on PreP, stop PreP
		if (pPerson->hiv().isOnPreP()) {
			pPerson->hiv().stopPreP();
		}

		// Schedule an initial monitoring event right away! (the 'true' is for 'right away')
		EventMonitoring *pEvtMonitor = new EventMonitoring(pPerson, true);
		population.onNewEvent(pEvtMonitor);
	} else {
		// If not on PreP, schedule new test.
		if (!pPerson->hiv().isOnPreP()) {
			EventHIVTest *pEvt = new EventHIVTest(pPerson);
			population.onNewEvent(pEvt);
		}
	}
}

double EventHIVTest::calculateInternalTimeInterval(const State *pState, double t0, double dt)
{
	// This is for the diagnosis event that should be scheduled right after the
	// screening event
	if (m_scheduleImmediately)
	{
		double hour = 1.0/(365.0*24.0); // an hour in a unit of a year
		return hour;
	}

	Person *pPerson = getPerson(0);
	double tMax = getTMax(pPerson);

	HazardFunctionHIVTest h0(pPerson, s_baseline, s_ageFactor, s_genderFactor, s_diagPartnersFactor, s_numPartnersFactor,
			           s_isDiagnosedFactor, s_healthSeekingPropensityFactor, s_beta, s_HSV2factor);
	TimeLimitedHazardFunction h(h0, tMax);

	return h.calculateInternalTimeInterval(t0, dt);
}

double EventHIVTest::solveForRealTimeInterval(const State *pState, double Tdiff, double t0)
{
	// This is for the diagnosis event that should be scheduled right after the
	// screening event
	if (m_scheduleImmediately)
	{
		double hour = 1.0/(365.0*24.0); // an hour in a unit of a year
		return hour;
	}

	Person *pPerson = getPerson(0);
	double tMax = getTMax(pPerson);

	HazardFunctionHIVTest h0(pPerson, s_baseline, s_ageFactor, s_genderFactor, s_diagPartnersFactor, s_numPartnersFactor,
			           s_isDiagnosedFactor, s_healthSeekingPropensityFactor, s_beta, s_HSV2factor);
	TimeLimitedHazardFunction h(h0, tMax);

	return h.solveForRealTimeInterval(t0, Tdiff);
}

double EventHIVTest::getTMax(const Person *pPerson)
{
	assert(pPerson != 0);
	double tb = pPerson->getDateOfBirth();

	assert(s_tMax > 0);
	return tb + s_tMax;
}

double EventHIVTest::s_baseline = 0;
double EventHIVTest::s_ageFactor = 0;
double EventHIVTest::s_genderFactor = 0;
double EventHIVTest::s_diagPartnersFactor = 0;
double EventHIVTest::s_numPartnersFactor = 0;
double EventHIVTest::s_isDiagnosedFactor = 0;
double EventHIVTest::s_healthSeekingPropensityFactor = 0;
double EventHIVTest::s_beta = 0;
double EventHIVTest::s_HSV2factor = 0;
double EventHIVTest::s_tMax = 0;


void EventHIVTest::processConfig(ConfigSettings &config, GslRandomNumberGenerator *pRndGen)
{
	bool_t r;

	if (!(r = config.getKeyValue("hivtest.baseline", s_baseline)) ||
	    !(r = config.getKeyValue("hivtest.agefactor", s_ageFactor)) ||
	    !(r = config.getKeyValue("hivtest.genderfactor", s_genderFactor)) ||
	    !(r = config.getKeyValue("hivtest.diagpartnersfactor", s_diagPartnersFactor)) ||
		!(r = config.getKeyValue("hivtest.numpartnersfactor", s_numPartnersFactor)) ||
	    !(r = config.getKeyValue("hivtest.isdiagnosedfactor", s_isDiagnosedFactor)) ||
		!(r = config.getKeyValue("hivtest.healthseekingpropensityfactor", s_healthSeekingPropensityFactor)) ||
	    !(r = config.getKeyValue("hivtest.beta", s_beta)) ||
	    !(r = config.getKeyValue("hivtest.t_max", s_tMax))||
	    !(r = config.getKeyValue("hivtest.HSV2factor", s_HSV2factor))
	   )
		abortWithMessage(r.getErrorString());
}

void EventHIVTest::obtainConfig(ConfigWriter &config)
{
	bool_t r;

	if (!(r = config.addKey("hivtest.baseline", s_baseline)) ||
	    !(r = config.addKey("hivtest.agefactor", s_ageFactor)) ||
	    !(r = config.addKey("hivtest.genderfactor", s_genderFactor)) ||
	    !(r = config.addKey("hivtest.diagpartnersfactor", s_diagPartnersFactor)) ||
		!(r = config.addKey("hivtest.numpartnersfactor", s_numPartnersFactor)) ||
	    !(r = config.addKey("hivtest.isdiagnosedfactor", s_isDiagnosedFactor)) ||
		!(r = config.addKey("hivtest.healthseekingpropensityfactor", s_healthSeekingPropensityFactor)) ||
	    !(r = config.addKey("hivtest.beta", s_beta)) ||
	    !(r = config.addKey("hivtest.t_max", s_tMax)) ||
	    !(r = config.addKey("hivtest.HSV2factor", s_HSV2factor))
	   )
		abortWithMessage(r.getErrorString());
}

// exp(baseline + ageFactor*(t-t_birth) + genderFactor*gender + diagPartnersFactor*numDiagnosedPartners + numPartnersFactor*numPartners +
//     isDiagnosedFactor*hasBeenDiagnosed + healthSeekingPropensityFactor*healthSeekingPropensity + HSV2factor*HSV2 + beta*(t-t_last_test))
//
// = exp(A + B*t) with
// 
//  A = baseline - ageFactor*t_birth + genderFactor*gender + diagPartnersFactor*numDiagnosedPartners
//  + numPartnersFactor*numRelationships + HSV2factor*HSV2
//      + isDiagnosedFactor*hasBeenDiagnosed + healthSeekingPropensityFactor*healthSeekingPropensity - beta*t_infected
//  B = ageFactor + beta
HazardFunctionHIVTest::HazardFunctionHIVTest(Person *pPerson, double baseline, double ageFactor,
		                                                 double genderFactor, double diagPartnersFactor, double numPartnersFactor,
							         double isDiagnosedFactor, double healthSeekingPropensityFactor, double beta, double HSV2factor)
	: m_baseline(baseline), m_ageFactor(ageFactor), m_genderFactor(genderFactor),
	  m_diagPartnersFactor(diagPartnersFactor), m_numPartnersFactor(numPartnersFactor), m_isDiagnosedFactor(isDiagnosedFactor), m_healthSeekingPropensityFactor(healthSeekingPropensityFactor),
	  m_beta(beta), m_HSV2factor(HSV2factor)
{
	assert(pPerson != 0);
	m_pPerson = pPerson;

	double tb = pPerson->getDateOfBirth();
	double ttest = pPerson->hiv().getLastTestTime();
	double G = (pPerson->isMan())?0:1;
	int D = pPerson->getNumberOfDiagnosedPartners();
	int P = pPerson->getNumberOfRelationships();
	int hasBeenDiagnosed = (pPerson->hiv().isDiagnosed())?1:0;
	double H = pPerson->getHealthSeekingPropensity();
	int HSV2 = (pPerson->hsv2().isInfected())?1:0;

	double A = baseline - ageFactor*tb + genderFactor*G + diagPartnersFactor*D + numPartnersFactor*P + isDiagnosedFactor*hasBeenDiagnosed + healthSeekingPropensityFactor*H - beta*ttest + HSV2factor*HSV2;
	double B = ageFactor + beta;

	setAB(A, B);
}

// This implementation is not necessary for running, it is provided for testing purposes
double HazardFunctionHIVTest::evaluate(double t)
{
	double tb = m_pPerson->getDateOfBirth();
	double ttest = m_pPerson->hiv().getLastTestTime();
	double G = (m_pPerson->isMan())?0:1;
	int D = m_pPerson->getNumberOfDiagnosedPartners();
	int P = m_pPerson->getNumberOfRelationships();
	int hasBeenDiagnosed = (m_pPerson->hiv().isDiagnosed())?1:0;
	double H = m_pPerson->getHealthSeekingPropensity();
	int HSV2 = (m_pPerson->hsv2().isInfected()) ?1:0;

	double age = (t-tb);

	return std::exp(m_baseline + m_ageFactor*age + m_genderFactor*G + m_diagPartnersFactor*D + m_numPartnersFactor*P +
			m_isDiagnosedFactor*hasBeenDiagnosed + m_healthSeekingPropensityFactor*H + m_beta*(t-ttest)+ m_HSV2factor*HSV2);
}

ConfigFunctions hivtestConfigFunctions(EventHIVTest::processConfig, EventHIVTest::obtainConfig, "EventHIVTest");

JSONConfig hivtestJSONConfig(R"JSON(
        "EventHIVTest": {
            "depends": null,
            "params": [
                [ "hivtest.baseline", 0 ],
                [ "hivtest.agefactor", 0 ],
                [ "hivtest.genderfactor", 0 ],
                [ "hivtest.diagpartnersfactor", 0 ],
				[ "hivtest.numpartnersfactor", 0 ],
                [ "hivtest.isdiagnosedfactor", 0 ],
                [ "hivtest.beta", 0 ],
				[ "hivtest.healthseekingpropensityfactor", 0 ],
		     	[ "hivtest.HSV2factor", 0 ],
                [ "hivtest.t_max", 200 ]	
            ],
            "info": [
                "When a person gets infected or drops out of treatment, a diagnosis event is ",
                "scheduled of which the fire time is determined by the following hazard:",
                "",
                " h = exp(baseline + agefactor*A(t) + genderfactor*G ",
                "         + diagpartnersfactor*ND + numpartnersfactor*P + isdiagnosedfactor*D",
                "         + beta*t + HSV2factor*HSV2)",
                "",
                "Here, A(t) is the age of the person, G is the gender (0 for a man, 1 for a",
                "woman), ND is the number of diagnosed partners, P is the number of partners and D is a flag (0 or 1)",
                "indicating if the person has been on treatment before (to have different",
                "behaviour for first diagnosis and re-testing after dropout)",
				"and H is the health-seeking propensity of the person."
            ]
        })JSON");


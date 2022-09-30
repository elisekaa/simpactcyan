#include "eventprepscreening.h"

#include "configdistributionhelper.h"
#include "configfunctions.h"
#include "configsettings.h"
#include "configwriter.h"
#include "eventhivtest.h"
#include "gslrandomnumbergenerator.h"
#include "jsonconfig.h"
#include "probabilitydistribution.h"

using namespace std;

EventPrePScreening::EventPrePScreening(Person *pPerson, bool scheduleImmediately): SimpactEvent(pPerson), m_scheduleImmediately(scheduleImmediately) {}

EventPrePScreening::~EventPrePScreening() {}

string EventPrePScreening::getDescription(double tNow) const
{
	return strprintf("PreP screening event for %s", getPerson(0)->getName().c_str());
}

void EventPrePScreening::writeLogs(const SimpactPopulation &pop, double tNow) const
{
	Person *pPerson = getPerson(0);
	writeEventLogStart(false, "prepscreening", tNow, pPerson, 0);
}

void EventPrePScreening::fire(Algorithm *pAlgorithm, State *pState, double t)
{
	SimpactPopulation &population = SIMPACTPOPULATION(pState);
	Person *pPerson = getPerson(0);

	// Immediately schedule test if not yet diagnosed
	if (!pPerson->hiv().isDiagnosed()) {
		EventHIVTest *pEvt = new EventHIVTest(pPerson, true);
		population.onNewEvent(pEvt);
	}

	// TODO if STI infected, immediately schedule treatment

	// If not HIV infected, schedule new screening
	if (!pPerson->hiv().isInfected()) {
		EventPrePScreening *pNewScreening = new EventPrePScreening(pPerson);
		population.onNewEvent(pNewScreening);
	}
}

void EventPrePScreening::processConfig(ConfigSettings &config, GslRandomNumberGenerator *pRndGen)
{
	delete s_pScreeningIntervalDistribution;
	s_pScreeningIntervalDistribution = getDistributionFromConfig(config, pRndGen, "prepscreening.interval");
}

void EventPrePScreening::obtainConfig(ConfigWriter &config)
{
	assert(s_pScreeningIntervalDistribution);
	addDistributionToConfig(s_pScreeningIntervalDistribution, config, "prepscreening.interval");
}

bool EventPrePScreening::isUseless(const PopulationStateInterface &pop)
{
	// PreP screening event becomes useless if person has been diagnosed with HIV
	Person *pPerson = getPerson(0);

	if (pPerson->hiv().isDiagnosed()) {
		return true;
	}

	return false;
}

double EventPrePScreening::getNewInternalTimeDifference(GslRandomNumberGenerator *pRndGen, const State *pState)
{
	// This is for the screening event that should be scheduled right after the
	// prep start event event
	if (m_scheduleImmediately)
	{
		double hour = 1.0/(365.0*24.0); // an hour in a unit of a year
		//return hour * pRndGen->pickRandomDouble();
		return hour;
	}

	double dt = s_pScreeningIntervalDistribution->pickNumber();

	assert(dt >= 0);

	return dt;
}

ProbabilityDistribution *EventPrePScreening::s_pScreeningIntervalDistribution = 0;

ConfigFunctions prepScreeningConfigFunctions(EventPrePScreening::processConfig, EventPrePScreening::obtainConfig, "EventPrePScreening");

JSONConfig prepScreeningJSONConfig(R"JSON(
        "EventPrePScreening": {
            "depends": null,
            "params": [ ["prepscreening.interval.dist", "distTypes", ["fixed", [ [ "value", 0.25 ] ] ] ] ],
            "info": [
                "TODO"
            ]
        })JSON");

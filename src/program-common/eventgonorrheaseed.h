#ifndef EVENTGONORRHEASEED_H
#define EVENTGONORRHEASEED_H

#include "eventseedbase.h"

class ConfigSettings;

class EventGonorrheaSeed : public EventSeedBase
{
public:
	EventGonorrheaSeed();
	~EventGonorrheaSeed();

	std::string getDescription(double tNow) const;
	void writeLogs(const SimpactPopulation &pop, double tNow) const;

	void fire(Algorithm *pAlgorithm, State *pState, double t);

	static void processConfig(ConfigSettings &config, GslRandomNumberGenerator *pRndGen);
	static void obtainConfig(ConfigWriter &config);

	static double getSeedTime() { return s_settings.m_seedTime; }

private:
	double getNewInternalTimeDifference(GslRandomNumberGenerator *pRndGen, const State *pState);

	static SeedEventSettings s_settings;
};

#endif // EVENTGONORRHEASEED_H

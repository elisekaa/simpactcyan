#include "person_gonorrhea.h"

Person_Gonorrhea::Person_Gonorrhea(Person *pSelf): Person_STI(pSelf), m_diseaseStage(Susceptible)
{}

Person_Gonorrhea::~Person_Gonorrhea() {}

bool Person_Gonorrhea::isInfected() const
{
	if (m_diseaseStage == Asymptomatic || m_diseaseStage == Symptomatic) {
		return true;
	} else {
		return false;
	}
}

void Person_Gonorrhea::setInfected(double t, Person *pOrigin, InfectionType iType)
{
	assert(iType != None);
	assert(!(pOrigin == 0 && iType != Seed));

	m_infectionTime = t;
	m_pInfectionOrigin = pOrigin;
	m_infectionType = iType;

	m_diseaseStage = Asymptomatic;
}

void Person_Gonorrhea::setRecovered(double t)

{
	m_diseaseStage = Susceptible;

	// Reset disease variables
	m_infectionTime = -1e200;
	m_pInfectionOrigin = 0;
	m_infectionType = None;
}

/*
 *
#include "person.h"
#include "configsettings.h"
#include "configwriter.h"
#include "configdistributionhelper.h"
#include "configfunctions.h"
#include "jsonconfig.h"
#include <vector>
#include <iostream>

using namespace std;

void Person_HSV2::processConfig(ConfigSettings &config, GslRandomNumberGenerator *pRndGen)
{
	assert(pRndGen != 0);

	delete m_pADist;
	delete m_pB2Dist;
	m_pADist = getDistributionFromConfig(config, pRndGen, "person.hsv2.a");
	m_pB2Dist = getDistributionFromConfig(config, pRndGen, "person.hsv2.b2");
}

void Person_HSV2::obtainConfig(ConfigWriter &config)
{
	addDistributionToConfig(m_pADist, config, "person.hsv2.a");
	addDistributionToConfig(m_pB2Dist, config, "person.hsv2.b2");
}

ConfigFunctions personHSVConfigFunctions(Person_HSV2::processConfig, Person_HSV2::obtainConfig, "Person_HSV2");

JSONConfig personHSV2JSONConfig(R"JSON(

        "PersonHSV2": {
            "depends": null,
            "params": [
				[ "person.hsv2.a.dist", "distTypes", [ "fixed", [ [ "value", 0 ]   ] ] ],
				[ "person.hsv2.b2.dist", "distTypes", [ "fixed", [ [ "value", 0 ]   ] ] ]
            ],
            "info": [
				"The 'a' parameter in the HSV2 transmission hazard is chosen from this",
				"distribution, allowing transmission to depend more on the individual",
				"The 'b2' parameter in the HSV2 transmission hazard is chosen from this",
				"distribution, allowing transmission to",
				"depend more on susceptibility for HSV2 only."
            ]
        })JSON");

 *
 */
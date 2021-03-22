library("irace")

# No utilizar más de 25 núcleos en la parametrización.
scenario <- readScenario(filename = "tuning/scenario.txt",
                         scenario = defaultScenario())

irace.main(scenario = scenario)
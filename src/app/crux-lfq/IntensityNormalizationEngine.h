#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "CQStatistics.h"
#include "Results.h"
#include "io/carp.h"

using std::vector;

namespace CruxQuant {

enum class Transform {
    Linear,
    Logarithmic,
    ExponentialAverage
};

// IParameterSampler interface
class IParameterSampler {
   public:
    virtual double Sample(double min, double max) = 0;
};

class ITransform {
   public:
    virtual double Transform(double min, double max, IParameterSampler& sampler) = 0;
};

class LinearTransform : public ITransform {
   public:
    double Transform(double min, double max, IParameterSampler& sampler) override {
        return sampler.Sample(min, max);
    }
};

class LogarithmicTransform : public ITransform {
   public:
    double Transform(double min, double max, IParameterSampler& sampler) override {
        if (min <= 0.0 || max <= 0.0) {
            throw std::invalid_argument("logarithmic scale requires min: " + std::to_string(min) + " and max: " + std::to_string(max) + " to be larger than zero");
        }

        double min2 = std::log10(min);
        double max2 = std::log10(max);
        double y = sampler.Sample(min2, max2);
        return std::pow(10.0, y);
    }
};

class ExponentialAverageTransform : public ITransform {
   public:
    double Transform(double min, double max, IParameterSampler& sampler) override {
        if (min >= 1.0 || max >= 1.0) {
            throw std::invalid_argument("ExponentialAverage scale requires min: " + std::to_string(min) + " and max: " + std::to_string(max) + " to be smaller than one");
        }

        double min2 = std::log10(1.0 - max);
        double max2 = std::log10(1.0 - min);
        double y = sampler.Sample(min2, max2);
        return 1.0 - std::pow(10.0, y);
    }
};

class TransformFactory {
   public:
    static ITransform* Create(Transform transform) {
        switch (transform) {
            case Transform::Linear:
                return new LinearTransform();
            case Transform::Logarithmic:
                return new LogarithmicTransform();
            case Transform::ExponentialAverage:
                return new ExponentialAverageTransform();
            default:
                throw std::invalid_argument("Unsupported transform: " + static_cast<int>(transform));
        }
    }
};

class ParameterBounds {
   public:
    double Min;
    double Max;

   private:
    ITransform* m_transform;

   public:
    ParameterBounds() : Min(0), Max(0), m_transform(TransformFactory::Create(Transform::Linear)) {}
    ParameterBounds(double min, double max, Transform transform = Transform::Linear)
        : Min(min), Max(max) {
        if (min >= max) {
            throw std::invalid_argument("min: " + std::to_string(min) + " is larger than or equal to max: " + std::to_string(max));
        }

        m_transform = TransformFactory::Create(transform);
    }

    ParameterBounds(double min, double max, ITransform* transform)
        : Min(min), Max(max), m_transform(transform) {
        if (min >= max) {
            throw std::invalid_argument("min: " + std::to_string(min) + " is larger than or equal to max: " + std::to_string(max));
        }

        if (transform == nullptr) {
            throw std::invalid_argument("transform");
        }
    }

    double NextValue(IParameterSampler& sampler) {
        return m_transform->Transform(Min, Max, sampler);
    }

    ParameterBounds& operator=(const ParameterBounds& other) {
        if (this != &other) {
            m_transform = other.m_transform;
            Min = other.Min;
            Max = other.Max;
        }
        return *this;
    }
};

class OptimizerResult {
   public:
    double Error;
    std::vector<double> ParameterSet;

    OptimizerResult(std::vector<double> parameterSet, double error)
        : Error(error), ParameterSet(std::move(parameterSet)) {
        if (ParameterSet.empty()) {
            throw std::invalid_argument("parameterSet");
        }
    }
    OptimizerResult& operator=(const OptimizerResult& other) {
        if (this != &other) {
            Error = other.Error;
            ParameterSet = other.ParameterSet;
        }
        return *this;
    }
};

class IOptimizer {
   public:
    virtual OptimizerResult OptimizeBest(std::function<OptimizerResult(std::vector<double>)> functionToMinimize) = 0;
    virtual std::vector<OptimizerResult> Optimize(std::function<OptimizerResult(std::vector<double>)> functionToMinimize) = 0;
};

class NelderMeadWithStartPoints : public IOptimizer {
    int m_maxIterationsPrRestart;
    int m_maxIterationsWithoutImprovement;
    int m_maxRestarts;
    double m_alpha;
    double m_gamma;
    double m_rho;
    double m_sigma;
    double m_noImprovementThreshold;
    std::vector<ParameterBounds> m_parameters;
    std::mt19937 m_random;
    int m_maxFunctionEvaluations;
    std::vector<double> startingValue;
    int m_totalFunctionEvaluations;

   public:
    NelderMeadWithStartPoints(std::vector<ParameterBounds> parameters, std::vector<double> startingValue, int maxRestarts = 8, double noImprovementThreshold = 0.001,
                              int maxIterationsWithoutImprovement = 5, int maxIterationsPrRestart = 0, int maxFunctionEvaluations = 0,
                              double alpha = 1, double gamma = 2, double rho = -0.5, double sigma = 0.5) {
        if (parameters.empty()) {
            throw std::invalid_argument("parameters");
        }
        if (maxIterationsWithoutImprovement <= 0) {
            throw std::invalid_argument("maxIterationsWithoutImprovement must be at least 1");
        }
        if (maxFunctionEvaluations < 0) {
            throw std::invalid_argument("maxFunctionEvaluations must be at least 1");
        }

        m_maxRestarts = maxRestarts;
        m_maxIterationsPrRestart = maxIterationsPrRestart;
        m_alpha = alpha;
        m_gamma = gamma;
        m_rho = rho;
        m_sigma = sigma;
        m_parameters = parameters;
        m_noImprovementThreshold = noImprovementThreshold;
        m_maxIterationsWithoutImprovement = maxIterationsWithoutImprovement;
        m_maxFunctionEvaluations = maxFunctionEvaluations;
        this->startingValue = startingValue;

        m_random = std::mt19937(startingValue[0]);  // Assuming startingValue[0] can be used as a seed
    }

    OptimizerResult OptimizeBest(std::function<OptimizerResult(std::vector<double>)> functionToMinimize) override {
        return Optimize(functionToMinimize).front();
    }

    std::vector<OptimizerResult> Optimize(std::function<OptimizerResult(std::vector<double>)> functionToMinimize) override {
        auto dim = m_parameters.size();
        std::vector<OptimizerResult> allResults;
        m_totalFunctionEvaluations = 0;
        vector<double> initialPoint(dim);
        for (int i = 0; i < initialPoint.size(); i++) {
            initialPoint[i] = startingValue[i];
        }
        for (int restarts = 0; restarts < m_maxRestarts; restarts++) {
            auto prevBest = EvaluateFunction(functionToMinimize, initialPoint);
            auto iterationsWithoutImprovement = 0;
            std::vector<OptimizerResult> results{prevBest};

            for (int i = 0; i < dim; i++) {
                auto a = (0.02 + 0.08 * m_random()) * (m_parameters[i].Max - m_parameters[i].Min);  // % simplex size between 2%-8% of min(xrange)

                auto p = a * (std::sqrt(dim + 1) + dim - 1) / (dim * std::sqrt(2));
                auto q = a * (std::sqrt(dim + 1) - 1) / (dim * std::sqrt(2));

                auto x = initialPoint;  // In C++, this will create a copy of initialPoint
                x[i] = x[i] + p;

                for (int j = 0; j < dim; j++) {
                    if (j != i) {
                        x[j] = x[j] + q;
                    }
                }

                BoundCheck(x);
                auto score = EvaluateFunction(functionToMinimize, x);
                results.push_back(score);
            }

            // simplex iter
            auto iterations = 0;
            while (true) {
                std::sort(
                    results.begin(),
                    results.end(),
                    [](const OptimizerResult& a, const OptimizerResult& b) {
                        return a.Error < b.Error;
                    });
                auto best = results.front();

                // break after m_maxIterationsPrRestart
                if (iterations >= m_maxIterationsPrRestart && m_maxIterationsPrRestart != 0) {
                    allResults.insert(allResults.end(), results.begin(), results.end());
                    break;
                }

                iterations++;
                double percentImprovement = -((best.Error - prevBest.Error) / prevBest.Error);

                if (percentImprovement > m_noImprovementThreshold) {
                    iterationsWithoutImprovement = 0;
                    prevBest = best;
                } else {
                    iterationsWithoutImprovement++;
                }

                // break after no_improv_break iterations with no improvement
                if (iterationsWithoutImprovement >= m_maxIterationsWithoutImprovement) {
                    allResults.insert(allResults.end(), results.begin(), results.end());
                    break;
                }

                // check if m_maxFunctionEvaluations is reached
                if (m_totalFunctionEvaluations >= m_maxFunctionEvaluations && m_maxFunctionEvaluations != 0) {
                    allResults.insert(allResults.end(), results.begin(), results.end());
                    break;
                }

                // centroid
                std::vector<double> x0(dim);

                for (auto it = results.begin(); it != results.end() - 1; ++it) {
                    auto parameters = it->ParameterSet;
                    for (int i = 0; i < parameters.size(); i++) {
                        x0[i] += parameters[i] / static_cast<double>(results.size() - 1);
                    }
                }

                BoundCheck(x0);

                // reflection
                auto last = results.back();
                std::vector<double> xr = x0;
                for (int i = 0; i < xr.size(); i++) {
                    xr[i] += (x0[i] - last.ParameterSet[i]) * m_alpha;
                }
                BoundCheck(xr);
                auto reflectionScore = EvaluateFunction(functionToMinimize, xr);

                double first = results.front().Error;
                if (first <= reflectionScore.Error && reflectionScore.Error < results[results.size() - 2].Error) {
                    results.pop_back();
                    results.push_back(reflectionScore);
                    continue;
                }

                // expansion
                if (reflectionScore.Error < first) {
                    std::vector<double> xe = x0;
                    for (int i = 0; i < xe.size(); i++) {
                        xe[i] += (x0[i] - last.ParameterSet[i]) * m_gamma;
                    }
                    BoundCheck(xe);
                    auto expansionScore = EvaluateFunction(functionToMinimize, xe);
                    if (expansionScore.Error < reflectionScore.Error) {
                        results.pop_back();
                        results.push_back(expansionScore);
                        continue;
                    } else {
                        results.pop_back();
                        results.push_back(reflectionScore);
                        continue;
                    }
                }

                // contraction
                std::vector<double> xc = x0;
                for (int i = 0; i < xc.size(); i++) {
                    xc[i] += (x0[i] - last.ParameterSet[i]) * m_rho;
                }
                BoundCheck(xc);
                auto contractionScore = EvaluateFunction(functionToMinimize, xc);
                if (contractionScore.Error < last.Error) {
                    results.pop_back();
                    results.push_back(contractionScore);
                    continue;
                }

                // shrink
                std::vector<double> x1 = results[0].ParameterSet;
                std::vector<OptimizerResult> nres;
                for (auto& tup : results) {
                    std::vector<double> xs = x1;
                    for (int i = 0; i < xs.size(); i++) {
                        xs[i] += (x1[i] - tup.ParameterSet[i]) * m_sigma;
                    }
                    BoundCheck(xs);
                    auto score = EvaluateFunction(functionToMinimize, xs);
                    nres.push_back(score);
                }

                results = nres;
            }
            // check if m_maxFunctionEvaluations is reached
            if (m_totalFunctionEvaluations >= m_maxFunctionEvaluations && m_maxFunctionEvaluations != 0) {
                allResults.insert(allResults.end(), results.begin(), results.end());
                break;
            }
        }
        std::vector<OptimizerResult> validResults;
        std::copy_if(
            allResults.begin(),
            allResults.end(),
            std::back_inserter(validResults),
            [](const OptimizerResult& r) {
                return !std::isnan(r.Error);
            });
        std::sort(
            validResults.begin(),
            validResults.end(),
            [](const OptimizerResult& a, const OptimizerResult& b) {
                return a.Error < b.Error;
            });
    }

   private:
    OptimizerResult EvaluateFunction(std::function<OptimizerResult(std::vector<double>)> functionToMinimize, std::vector<double> parameters) {
        m_totalFunctionEvaluations++;
        return functionToMinimize(parameters);
    }

    void BoundCheck(std::vector<double>& parameters) {
        for (int i = 0; i < parameters.size(); i++) {
            auto parameter = m_parameters[i];
            parameters[i] = std::max(parameter.Min, std::min(parameters[i], parameter.Max));
        }
    }
};

class IntensityNormalizationEngine {
   private:
    static const int numPeptidesDesiredFromEachFraction = 500;
    static const int numPeptidesDesiredInMatrix = 5000;
    CruxLFQResults results;
    bool integrate;
    bool quantifyAmbiguousPeptides;

    void NormalizeTechreps() {
        std::vector<CruxQuant::Peptides> peptides;
        peptides.reserve(results.PeptideModifiedSequences.size());

        std::transform(
            results.PeptideModifiedSequences.begin(),
            results.PeptideModifiedSequences.end(),
            std::back_inserter(peptides),
            [](const std::pair<const std::string, CruxQuant::Peptides>& pair) {
                return pair.second;
            });

        std::unordered_map<std::string, std::vector<SpectraFileInfo>> conditions;

        for (const auto& spectraFile : results.spectraFiles) {
            conditions[spectraFile.Condition].push_back(spectraFile);
        }

        for (const auto& condition : conditions) {
            std::unordered_map<int, std::vector<SpectraFileInfo>> bioreps;

            for (const auto& spectraFile : condition.second) {
                bioreps[spectraFile.BiologicalReplicate].push_back(spectraFile);
            }

            for (const auto& biorep : bioreps) {
                std::unordered_map<int, std::vector<SpectraFileInfo>> fractions;

                for (const auto& spectraFile : biorep.second) {
                    fractions[spectraFile.Fraction].push_back(spectraFile);
                }

                for (const auto& fraction : fractions) {
                    const auto& techreps = fraction.second;

                    for (size_t t = 1; t < techreps.size(); ++t) {
                        std::vector<double> foldChanges;

                        for (size_t p = 0; p < peptides.size(); ++p) {
                            double techrep1Intensity = peptides[p].getIntensity(techreps[0].FullFilePathWithExtension);
                            double techrepTIntensity = peptides[p].getIntensity(techreps[t].FullFilePathWithExtension);

                            if (techrep1Intensity > 0 && techrepTIntensity > 0) {
                                foldChanges.push_back(techrepTIntensity / techrep1Intensity);
                            }
                        }

                        if (foldChanges.empty()) {
                            // TODO: throw an exception?
                            return;
                        }

                        double medianFoldChange = Median(foldChanges);
                        double normalizationFactor = 1.0 / medianFoldChange;

                        // normalize to median fold-change
                        for (auto& peak : results.Peaks[techreps[t].FullFilePathWithExtension]) {
                            for (auto& isotopeEnvelope : peak.isotopicEnvelopes) {
                                isotopeEnvelope.Normalize(normalizationFactor);
                            }

                            // recalculate intensity after normalization
                            peak.calculateIntensityForThisFeature(integrate);
                        }
                    }
                }
            }
        }
    }

    void NormalizeFractions() {
        auto& spectraFiles = results.spectraFiles;
        auto& peaks = results.Peaks;
        if (std::max_element(
                results.spectraFiles.begin(),
                results.spectraFiles.end(),
                [](const SpectraFileInfo& a, const SpectraFileInfo& b) {
                    return a.Fraction < b.Fraction;
                })
                ->Fraction == 0) {
            return;
        }
        vector<CruxQuant::Peptides> peptides;
        peptides.reserve(results.PeptideModifiedSequences.size());
        std::transform(
            results.PeptideModifiedSequences.begin(),
            results.PeptideModifiedSequences.end(),
            std::back_inserter(peptides),
            [](const std::pair<const std::string, CruxQuant::Peptides>& pair) {
                return pair.second;
            });

        std::set<std::string> conditionsSet;

        for (const auto& spectraFile : results.spectraFiles) {
            conditionsSet.insert(spectraFile.Condition);
        }

        std::vector<std::string> conditions(conditionsSet.begin(), conditionsSet.end());

        std::vector<SpectraFileInfo> filesForCond1Biorep1;
        for (const auto& spectraFile : results.spectraFiles) {
            if (spectraFile.Condition == conditions[0] &&
                spectraFile.BiologicalReplicate == 0 &&
                spectraFile.TechnicalReplicate == 0) {
                filesForCond1Biorep1.push_back(spectraFile);
            }
        }

        for (const auto& condition : conditions) {
            auto filesForThisCondition = std::vector<SpectraFileInfo>();
            std::copy_if(spectraFiles.begin(), spectraFiles.end(), std::back_inserter(filesForThisCondition),
                         [&condition](const SpectraFileInfo& p) { return p.Condition == condition; });
            auto numB = std::count_if(filesForThisCondition.begin(), filesForThisCondition.end(),
                                      [&condition](const SpectraFileInfo& p) { return p.BiologicalReplicate == 0; });

            for (int b = 0; b < numB; ++b) {
                if (b == 0 && condition == conditions[0]) {
                    continue;  // condition 1 biorep 1 is the reference, don't normalize it
                }

                carp(CARP_INFO, "Normalizing condition \"%s\" biorep %d", condition.c_str(), b + 1);

                auto filesForThisBiorep = std::vector<SpectraFileInfo>();
                std::copy_if(
                    filesForThisCondition.begin(),
                    filesForThisCondition.end(),
                    std::back_inserter(filesForThisBiorep),
                    [&b](const SpectraFileInfo& p) { return p.BiologicalReplicate == b && p.TechnicalReplicate == 0; });
                auto numF = std::max_element(
                                filesForCond1Biorep1.begin(),
                                filesForCond1Biorep1.end(),
                                [](const SpectraFileInfo& a, const SpectraFileInfo& b) { return a.Fraction < b.Fraction; })
                                ->Fraction +
                            1;

                auto seenInBothBioreps = std::vector<Peptides>();
                for (auto& peptide : peptides) {
                    bool seenInBiorep1 = false;
                    bool seenInBiorep2 = false;

                    for (const auto& file : filesForCond1Biorep1) {
                        if (peptide.getIntensity(file.FullFilePathWithExtension) > 0) {
                            seenInBiorep1 = true;
                        }
                    }

                    for (const auto& file : filesForThisBiorep) {
                        if (peptide.getIntensity(file.FullFilePathWithExtension) > 0) {
                            seenInBiorep2 = true;
                        }
                    }

                    if (seenInBiorep1 && seenInBiorep2) {
                        seenInBothBioreps.push_back(peptide);
                    }
                }

                auto numP = seenInBothBioreps.size();
                auto myIntensityArray = std::vector<std::vector<std::vector<double>>>(numP, std::vector<std::vector<double>>(2, std::vector<double>(numF)));

                for (std::size_t p = 0; p < numP; ++p) {
                    auto& peptide = seenInBothBioreps[p];

                    for (const auto& file : filesForCond1Biorep1) {
                        myIntensityArray[p][0][file.Fraction] = peptide.getIntensity(file.FullFilePathWithExtension);
                    }

                    for (const auto& file : filesForThisBiorep) {
                        myIntensityArray[p][1][file.Fraction] = peptide.getIntensity(file.FullFilePathWithExtension);
                    }
                }

                auto normFactors = GetNormalizationFactors(myIntensityArray, numP, numF);
                carp(CARP_INFO, "Normalization factors for condition \"%s\" biorep %d:", condition.c_str(), b + 1);

                for (const auto& spectraFile : filesForThisBiorep) {
                    for (auto& peak : peaks[spectraFile.FullFilePathWithExtension]) {
                        for (auto& isotopeEnvelope : peak.isotopicEnvelopes) {
                            isotopeEnvelope.Normalize(normFactors[spectraFile.Fraction]);
                        }

                        peak.calculateIntensityForThisFeature(integrate);
                    }
                }
            }
        }
    }

    void NormalizeBioreps() {
        vector<Peptides> peptides;
        for (const auto& v : results.PeptideModifiedSequences) {
            peptides.push_back(v.second);
        }

        std::vector<std::vector<SpectraFileInfo>> conditions;
        std::map<std::string, std::vector<SpectraFileInfo>> grouped;

        for (const auto& v : results.spectraFiles) {
            grouped[v.Condition].push_back(v);
        }

        for (const auto& group : grouped) {
            conditions.push_back(group.second);
        }

        vector<vector<double>> biorepIntensityPair(peptides.size(), vector<double>(2));

        std::vector<SpectraFileInfo> firstConditionFirstBiorep;
        for (const auto& v : conditions[0]) {
            if (v.BiologicalReplicate == 0 && v.TechnicalReplicate == 0) {
                firstConditionFirstBiorep.push_back(v);
            }
        }
        for (const auto& file : firstConditionFirstBiorep) {
            for (size_t p = 0; p < peptides.size(); ++p) {
                biorepIntensityPair[p][0] += peptides[p].getIntensity(file.FullFilePathWithExtension);
            }
        }
        for (auto& condition : conditions) {
            std::map<int, std::vector<SpectraFileInfo>> bioreps;
            for (auto& v : condition) {
                bioreps[v.BiologicalReplicate].push_back(v);
                for (auto& biorep : bioreps) {
                    for (size_t p = 0; p < peptides.size(); ++p) {
                        biorepIntensityPair[p][1] = 0;
                    }

                    std::map<int, std::vector<SpectraFileInfo>> fractions;
                    for (auto& v : biorep.second) {
                        fractions[v.Fraction].push_back(v);
                    }

                    for (auto& fraction : fractions) {
                        SpectraFileInfo firstTechrep;
                        for (auto& v : fraction.second) {
                            if (v.TechnicalReplicate == 0) {
                                firstTechrep = v;
                                break;
                            }
                        }

                        for (size_t p = 0; p < peptides.size(); ++p) {
                            biorepIntensityPair[p][1] += peptides[p].getIntensity(firstTechrep.FullFilePathWithExtension);
                        }
                    }

                    std::vector<double> foldChanges;

                    for (size_t p = 0; p < peptides.size(); ++p) {
                        if (biorepIntensityPair[p][0] > 0 && biorepIntensityPair[p][1] > 0) {
                            foldChanges.push_back(biorepIntensityPair[p][1] / biorepIntensityPair[p][0]);
                        }
                    }

                    if (foldChanges.empty()) {
                        // TODO: throw an exception?
                        return;
                    }

                    double medianFoldChange = Median(foldChanges);
                    double normalizationFactor = 1.0 / medianFoldChange;

                    // normalize to median fold-change
                    for (auto& file : biorep.second) {
                        for (auto& peak : results.Peaks[file.FullFilePathWithExtension]) {
                            for (auto& isotopeEnvelope : peak.isotopicEnvelopes) {
                                isotopeEnvelope.Normalize(normalizationFactor);
                            }

                            // recalculate intensity after normalization
                            peak.calculateIntensityForThisFeature(integrate);
                        }
                    }
                }
            }
        }
    }

    double CalculateNormalizationFactorError(vector<double>& reference, vector<vector<double>>& sampleToNormalize, vector<double>& normalizationFactors, int numP, int numF) {
        double totalError = 0;

        for (int p = 0; p < numP; ++p) {
            // Sum the intensities with the current normalization factors
            double normalizedReplicateIntensity = 0;
            for (int f = 0; f < numF; ++f) {
                normalizedReplicateIntensity += sampleToNormalize[p][f] * normalizationFactors[f];
            }

            // Calculate the peptide error
            double peptideError = log(normalizedReplicateIntensity) - log(reference[p]);

            totalError += peptideError;
        }

        return totalError;
    }

    vector<double> GetNormalizationFactors(std::vector<std::vector<std::vector<double>>> peptideIntensities, int numP, int numF) {
        std::vector<double> referenceSample(numP, 0.0);
        std::vector<std::vector<double>> sampleToNormalize(numP, std::vector<double>(numF, 0.0));
        for (int p = 0; p < numP; ++p) {
            for (int f = 0; f < numF; ++f) {
                referenceSample[p] += peptideIntensities[p][0][f];
                sampleToNormalize[p][f] = peptideIntensities[p][1][f];
            }
        }
        std::vector<double> bestNormFactors(numF, 1.0);

        std::vector<double> initialErrors(numP, 0.0);
        double bestError = CalculateNormalizationFactorError(referenceSample, sampleToNormalize, bestNormFactors, numP, numF);
        // constraint (normalization factors must be >0.3 and <3)
        std::vector<ParameterBounds> parameterArray(numF);
        for (int f = 0; f < numF; f++) {
            parameterArray[f] = ParameterBounds(0.3, 3, Transform::Linear);
        }
        // find approximate best starting area for each fraction normalization factor
        for (int f = 0; f < numF; f++) {
            double bestFractionError = std::numeric_limits<double>::infinity();
            double start = parameterArray[0].Min;
            double end = parameterArray[0].Max;
            std::vector<double> factors(numF, 1.0);

            for (double n = start; n <= end; n += 0.01) {
                factors[f] = round(n * 100) / 100;

                double error = CalculateNormalizationFactorError(referenceSample, sampleToNormalize, factors, numP, numF);

                if (error < bestFractionError) {
                    bestFractionError = error;
                    bestNormFactors[f] = factors[f];
                }
            }
        }

        // find the best normalization factors (minimize error)
        std::vector<double> errors(numP);

        // define minimization metric
        auto minimize = [&](std::vector<double> v) -> OptimizerResult {
            // calculate error with these normalization factors
            double candidateError = CalculateNormalizationFactorError(referenceSample, sampleToNormalize, v, numP, numF);

            return OptimizerResult(v, candidateError);
        };

        // create optimizer
        NelderMeadWithStartPoints optimizer(parameterArray, bestNormFactors, 10);
        OptimizerResult result = optimizer.OptimizeBest(minimize);

        double sampleError = result.Error;
        vector<double> normalizationFactors = result.ParameterSet;

        if (sampleError < bestError) {
            if (sampleError < bestError) {
                bestError = sampleError;
                bestNormFactors = normalizationFactors;
            }
        }

        return bestNormFactors;
    }

   public:
    IntensityNormalizationEngine(CruxLFQResults results, bool integrate, bool quantifyAmbiguousPeptides = false)
        : results(results), integrate(integrate), quantifyAmbiguousPeptides(quantifyAmbiguousPeptides) {}

    void NormalizeResults() {
        results.calculatePeptideResults(quantifyAmbiguousPeptides);

        carp(CARP_INFO, "Normalizing fractions");
        NormalizeFractions();
        results.calculatePeptideResults(quantifyAmbiguousPeptides);

        carp(CARP_INFO, "Normalizing bioreps and conditions");
        NormalizeBioreps();
        results.calculatePeptideResults(quantifyAmbiguousPeptides);

        carp(CARP_INFO, "Normalizing techreps");
        NormalizeTechreps();
        results.calculatePeptideResults(quantifyAmbiguousPeptides);
    }
};
}  // namespace CruxQuant
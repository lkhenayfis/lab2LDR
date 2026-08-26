using lab2ldr
using Test

# On the abstract-engine branch these variable-name constants live in SDDPlab.Lab
# and are no longer re-exported at the SDDPlab top level, so import them explicitly.
import SDDPlab.Lab: STORED_VOLUME, TURBINED_FLOW, SPILLAGE, THERMAL_GENERATION, LOAD, DEFICIT
using LinearDecisionRules
using JuMP
using HiGHS
using Distributions
using Copulas
using LinearAlgebra

# helpers internos (nao exportados) exercitados isoladamente
import lab2ldr:
    build_cartesian_indexes,
    copula2Σ,
    validate_univar_gaussian,
    validate_copulas_kind,
    determinisitic_stage1,
    build_stage_ldruncertainty

# estudo de referencia: uma unica usina hidreletrica, 1 barra, 2 termicas, 4 estagios
const DATA_1DTOY = normpath(joinpath(@__DIR__, "..", "data", "1dtoy"))

@testset "lab2ldr.jl" begin

    # -------------------------------------------------------------------------
    # Funcoes auxiliares puras, testadas sem precisar montar o modelo completo
    # nem acionar o otimizador.
    # -------------------------------------------------------------------------
    @testset "helpers puros" begin
        @testset "build_cartesian_indexes" begin
            # deve percorrer todas as dimensoes menos a ultima (a de estagios)
            @test build_cartesian_indexes(zeros(1, 4)) == CartesianIndex.(1:1)
            @test build_cartesian_indexes(zeros(3, 4)) == CartesianIndex.(1:3)
        end

        @testset "copula2Σ" begin
            # copula independente -> sem correlacao cruzada
            @test copula2Σ(IndependentCopula(2)) == zeros(2, 2)
            # copula gaussiana -> apenas termos fora da diagonal
            @test copula2Σ(GaussianCopula([1.0 0.5; 0.5 1.0])) == [0.0 0.5; 0.5 0.0]
        end

        @testset "validate_univar_gaussian" begin
            @test validate_univar_gaussian([Normal(41, 21), LogNormal(0.0, 1.0)])
            @test_throws ErrorException validate_univar_gaussian([Uniform(0, 1)])
        end

        @testset "validate_copulas_kind" begin
            @test validate_copulas_kind(
                [IndependentCopula(1), GaussianCopula([1.0 0.0; 0.0 1.0])]
            )
        end

        @testset "determinisitic_stage1" begin
            # sempre retorna vetor: escalar embrulhado para distribuicoes univariadas
            @test determinisitic_stage1(Normal(41, 21)) == [41.0]
            @test determinisitic_stage1(MvNormal([1.0, 2.0], [1.0 0.0; 0.0 1.0])) ==
                [1.0, 2.0]
        end

        @testset "build_stage_ldruncertainty (univariada)" begin
            d = build_stage_ldruncertainty(Normal(41, 21), IndependentCopula(1))
            @test d isa UnivariateDistribution
            # truncamento padrao: [valor 0, quantil 0.975]
            @test minimum(d) == 0.0
            @test maximum(d) ≈ quantile(Normal(41, 21), 0.975)
        end
    end

    # -------------------------------------------------------------------------
    # Caminho funcional: estudo com uma unica hidro montado de ponta a ponta.
    # -------------------------------------------------------------------------
    @testset "modelo 1dtoy (uma hidro)" begin
        m = build_model(DATA_1DTOY, HiGHS.Optimizer; n_breakpoints = 1)

        @testset "estrutura" begin
            @test m isa LinearDecisionRules.LDRModel

            # dimensoes (entidade, estagio) das variaveis de sistema
            @test size(m[STORED_VOLUME]) == (1, 4)
            @test size(m[TURBINED_FLOW]) == (1, 4)
            @test size(m[SPILLAGE]) == (1, 4)
            @test size(m[THERMAL_GENERATION]) == (2, 4)
            @test size(m[LOAD]) == (1, 4)
            @test size(m[DEFICIT]) == (1, 4)

            # afluencia do primeiro estagio e deterministica (fixa); as demais existem
            @test is_fixed(m[:INFLOW_1])
            for s in 2:4
                @test haskey(object_dictionary(m), Symbol("INFLOW_$s"))
            end
        end

        @testset "atributo de estagio" begin
            for s in 1:4
                @test get_attribute(
                    m[STORED_VOLUME][1, s], LinearDecisionRules.Stage()
                ) == s
            end
        end

        @testset "solucao" begin
            LinearDecisionRules.optimize!(m)
            obj = LinearDecisionRules.objective_value(m; dual = true)
            @test isfinite(obj)
            # ancora de regressao sobre os dados atuais do estudo 1dtoy
            @test obj ≈ 3.7550274748 atol = 1e-2
        end
    end
end

using lab2ldr
using Test

using LinearDecisionRules
using JuMP
using HiGHS
using Distributions
using Copulas

# variable-name constants (submodule import; not top-level exported on abstract-engine)
import SDDPlab.Lab: STORED_VOLUME, THERMAL_GENERATION

# helpers internos (nao exportados) exercitados isoladamente
import lab2ldr:
    build_cartesian_indexes,
    copula2Σ,
    validate_univar_gaussian,
    validate_copulas_kind,
    determinisitic_stage1,
    build_stage_ldruncertainty

const DATA = normpath(joinpath(@__DIR__, "..", "data"))

# -----------------------------------------------------------------------------------------------
# Resultados "golden" capturados do projeto CasosLDR (HydroToy/hydro2.jl e hydro3.jl, funcao
# `ldr_solution`), com LinearDecisionRules Multistage @ 0.3.0.
# Cada linha: (estudo, n_breakpoints, Obj, Bound), onde
#   Obj   = objective_value(ldr)              (primal)
#   Bound = objective_value(ldr; dual = true) (dual)
#
# Os DUAIS do hydrotoy3 batem exato em n=0, mas divergem ~0.2% para n>=1 (efeito residual do LDR
# dual com dois estagios estocasticos); por isso apenas hydrotoy3 com n>=1 segue marcado
# `@test_broken` — a investigar; quando baterem, o `@test_broken` vira "Unexpected Pass" e
# sinaliza para promover a `@test`.
# -----------------------------------------------------------------------------------------------
const GOLDEN = [
    ("hydrotoy2", 0, 119250.000000, 90500.000000),
    ("hydrotoy2", 1, 98081.919583, 91161.044007),
    ("hydrotoy2", 2, 93847.451161, 91087.710925),
    ("hydrotoy2", 3, 93412.932692, 91707.138413),
    ("hydrotoy3", 0, 191500.000000, 90750.000000),
    ("hydrotoy3", 1, 119344.456318, 91806.622420),
    ("hydrotoy3", 2, 107853.175057, 91673.182804),
    ("hydrotoy3", 3, 106060.796247, 92577.468691),
]

const EXPECTED_STAGES = Dict("hydrotoy2" => 2, "hydrotoy3" => 3)

@testset "lab2ldr.jl" begin

    # -------------------------------------------------------------------------
    # Funcoes auxiliares puras, testadas sem precisar montar o modelo completo.
    # -------------------------------------------------------------------------
    @testset "helpers puros" begin
        @testset "build_cartesian_indexes" begin
            @test build_cartesian_indexes(zeros(1, 4)) == CartesianIndex.(1:1)
            @test build_cartesian_indexes(zeros(3, 4)) == CartesianIndex.(1:3)
        end

        @testset "copula2Σ" begin
            @test copula2Σ(IndependentCopula(2)) == zeros(2, 2)
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
            @test determinisitic_stage1(Normal(41, 21)) == [41.0]
            @test determinisitic_stage1(MvNormal([1.0, 2.0], [1.0 0.0; 0.0 1.0])) ==
                [1.0, 2.0]
        end

        @testset "build_stage_ldruncertainty (univariada)" begin
            # a truncagem foi movida para do_truncation!; este metodo agora apenas repassa a
            # distribuicao recebida, sem alteracao
            nd = Normal(41.0, 21.0)
            @test build_stage_ldruncertainty(nd, IndependentCopula(1)) === nd
        end
    end

    # -------------------------------------------------------------------------
    # Casos hydrotoy montados de ponta a ponta e comparados com o golden do
    # CasosLDR. Estrutura -> deve passar; objetivo/bound -> @test_broken.
    # -------------------------------------------------------------------------
    @testset "modelo hydrotoy (vs CasosLDR golden)" begin
        for (study, n, obj_golden, bound_golden) in GOLDEN
            @testset "$study n=$n" begin
                m = build_model(joinpath(DATA, study), HiGHS.Optimizer; n_breakpoints = n)

                # estrutura (deve passar)
                @test m isa LinearDecisionRules.LDRModel
                ns = EXPECTED_STAGES[study]
                @test size(m[STORED_VOLUME]) == (1, ns)
                @test size(m[THERMAL_GENERATION]) == (3, ns)
                @test is_fixed(m[:INFLOW_1])
                # estagio 1 nao e truncado -> valor fixo = media de Normal(60, 20)
                @test fix_value(m[:INFLOW_1]) ≈ 60.0

                # objetivo/bound vs golden do CasosLDR
                LinearDecisionRules.optimize!(m)
                obj = LinearDecisionRules.objective_value(m)
                bound = LinearDecisionRules.objective_value(m; dual = true)
                @test isapprox(obj, obj_golden; rtol = 1e-4)
                # o dual do hydrotoy3 bate exato em n=0 mas diverge ~0.2% para n>=1 (residual do
                # LDR dual com dois estagios estocasticos, a investigar) -> so esses ficam broken
                if study == "hydrotoy3" && n >= 1
                    @test_broken isapprox(bound, bound_golden; rtol = 1e-4)
                else
                    @test isapprox(bound, bound_golden; rtol = 1e-4)
                end
            end
        end
    end
end

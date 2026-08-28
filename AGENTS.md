# AGENTS.md — working in the matRad repository

Guidance for AI coding agents and human contributors. It complements
[CONTRIBUTING.md](CONTRIBUTING.md) (bug reports, PR process, programming style and the
contribution checklist) and
[test/README.md](test/README.md) (test framework). When in doubt, those documents and the
existing code win over assumptions.

matRad is an open-source MATLAB/Octave toolkit for radiotherapy treatment planning
(photons, protons, carbon/helium ions, brachytherapy, VMAT). Development happens on the
`dev` branch; `master` holds releases. Pull requests target `dev`.

## Repository layout

| Path | Purpose |
| --- | --- |
| `matRad/` | **The source code.** Everything that is part of the toolkit lives here (`doseCalc/`, `optimization/`, `steering/`, `sequencing/`, `bioModels/`, `dicom/`, `IO/`, `gui/`, `util/`, ...). `MatRad_Config.m` is the singleton configuration/logging class. |
| `matRad.m`, `matRadGUI.m`, `matRad_rc.m` | Top-level entry points. `matRad_rc` initializes paths and returns the `MatRad_Config` instance; call it before anything else. |
| `test/` | MOxUnit test suite, mirroring the `matRad/` folder structure. `helper_*.m` files are shared test helpers, `testData/` holds fixtures. |
| `examples/` | Documented workflow examples (`matRad_exampleNN_*.m`). Most of them are executed as tests in CI (`test/autoExampleTest`). |
| `docs/` | Sphinx documentation (reStructuredText, `sphinxcontrib-matlabdomain` for the API reference). |
| `tools/` | Helper scripts not part of the toolkit proper (see `tools/README.md`). |
| `standalone/`, `matRad_buildStandalone.m` | Assets and build script for the compiled standalone application (MATLAB Compiler). |
| `submodules/` | Git submodules with the *sources* of external tools (MCsquare, ompMC, MOxUnit, MOcov, matlab2tikz, latexTable). Run `git submodule update --init` after cloning. |
| `thirdParty/` | Prebuilt **third-party software** and binaries shipped with matRad (IPOPT mex files per platform/Octave version, MCsquare and ompMC binaries and data). Do not edit third-party code in place; changes go upstream or into the corresponding submodule, and rebuilt binaries are placed here. Respect the individual licenses in these folders. |
| `matRad/phantoms/`, `matRad/basedata/`, `matRad/hluts/` | Bundled patient phantoms, machine base data and HU lookup tables. |
| `userdata/` | User-specific data, not for tracked code. |
| `.github/` | CI workflows (tests on MATLAB R2022b, latest MATLAB and Octave; docs build; coverage; packaging). |

## Architecture: procedural top-level API on top of an object-oriented core

matRad deliberately has two layers:

1. **Top-level procedural API** — the `matRad_*` functions directly in `matRad/`
   (`matRad_generateStf`, `matRad_calcDoseInfluence`, `matRad_calcDoseForward`,
   `matRad_fluenceOptimization`, `matRad_sequencing`, `matRad_directApertureOptimization`,
   `matRad_planAnalysis`, ...). They operate on the plain data structures `ct`, `cst`,
   `pln`, `stf`, `dij`, `resultGUI` and are what users and the examples call. Their job is
   to **instantiate and configure the correct algorithm objects from `pln`** and hand
   results back as structs; they should stay thin and hold no algorithmic logic themselves.
2. **Object-oriented core** — the algorithm class hierarchies in the subfolders:
   stf generators (`matRad_StfGeneratorBase`, configured from `pln.propStf`), dose engines
   (`DoseEngines.matRad_DoseEngineBase`, `pln.propDoseCalc`), sequencers (`matRad_SequencerBase`, `pln.propSeq`), scenario models (`matRad_ScenarioModel`, `pln.multScen`), biological models (`matRad_BiologicalModel`, `pln.bioModel`), optimization functions/objectives, ...

The wiring follows a common pattern: each base class offers a
`get<Thing>FromPln(pln)` (or `create`/`validate`) static factory that selects the concrete
class by its constant `shortName` (e.g. `pln.propDoseCalc.engine = 'HongPB'`), falls back to
`matRad_cfg.defaults.prop*` when nothing valid is given, and then copies the remaining
`pln.prop*` fields onto the object's properties with `matRad_assignPropertiesFromStruct`
(unknown fields produce a warning). Subclasses are **discovered automatically**
(`matRad_findSubclasses`), so adding a new engine/generator/sequencer usually means adding
one class with `shortName`, `possibleRadiationModes`, `isAvailable` and the abstract
methods implemented, plus (if it should be a default) an entry in `MatRad_Config.defaults`.
Users may also bypass `pln` entirely and pass a configured object, so keep both routes
working and consistent (mismatches are warned about, not silently overridden).

When contributing:

- Put logic into the class hierarchy; expose it through `pln.prop*` options and, if it is
  a new workflow step, through a thin top-level function.
- Every configurable property needs a sensible default (in the class and/or
  `MatRad_Config.defaults`) so that a minimal `pln` keeps working.
- Document new `pln.prop*` fields in `docs/` and in the class docstring.

## `MatRad_Config` — defaults and logging

`MatRad_Config` is a singleton; obtain it with
`matRad_cfg = MatRad_Config.instance();` at the top of any function/method that needs it
(`matRad_rc` returns the same instance). It provides:

- **Logging:** `dispInfo`, `dispWarning`, `dispError`, `dispDebug`,
  `dispDeprecationWarning` (all `fprintf`-style). Use these instead of
  `disp`/`fprintf`/`warning`/`error`; output is filtered by `matRad_cfg.logLevel`
  (1 errors ... 5 debug) and redirected consistently in the GUI, standalone and tests.
  `dispError` throws, so no code after it is reached.
- **Defaults:** `matRad_cfg.defaults` holds the default machines, biological models and all
  `propStf`/`propDoseCalc`/`propOpt`/`propSeq`/... parameters. Read defaults from there
  (`matRad_getFieldOrDefault` helps) rather than hard-coding values; add new defaults in
  `setDefaultProperties`. `setDefaultPropertiesForTesting` provides the reduced settings
  used by the test runner.
- **Environment:** `isOctave`/`isMatlab`, `matRadRoot`, `matRadSrcRoot`, version and
  GPU/parallel flags.

Do not store per-call state in the config object, and never call `MatRad_Config` methods
that reset it (`reset`, `setDefaultProperties*`) from library code.

## Environments and compatibility

- **Stable MATLAB target: R2022b.** Code must run on R2022b; do not rely on language or
  toolbox features introduced later. CI additionally runs the latest MATLAB release.
- **Octave compatibility whenever possible.** Octave 8 (8.4.0) is tested in CI. Prefer
  constructs that work in both environments. Use `matRad_cfg.isOctave` / `isMatlab`
  (from `MatRad_Config.instance()`) for the rare environment-specific branches.
- **Octave compatibility layer:** `matRad/util/octaveCompat/` contains `*Compat`
  wrappers (`matRad_roundCompat`, `matRad_gatherCompat`, `matRad_getPropsCompat`,
  `matRad_ispropCompat`, `matRad_underlyingTypeCompat`, ...). Use them instead of MATLAB-only
  call signatures, and add a new wrapper there (with a test) when you hit another incompatibility.
  Known pitfalls: two-argument `round(x,n)`, `arguments` blocks, `string` arrays in many
  places, `isprop` on class objects, some OOP features (e.g. property validation syntax).
- Toolbox dependencies must be optional and checked at runtime (e.g.
  `matRad_checkEnvImageProcessingRequirements`), with a clear error or fallback when missing.
- Never use relative paths to repository files. Build absolute paths at runtime
  (`matRad_cfg.matRadRoot`, `mfilename('fullpath')`, `which`) and, if you have to `cd`, restore
  the previous directory afterwards.

## Coding conventions

- **Object-oriented design is preferred for new contributions.** Follow the existing
  patterns: abstract base classes plus concrete implementations that are discovered
  automatically (dose engines `matRad_DoseEngineBase`, stf generators, sequencers
  `matRad_SequencerBase`, optimizers, scenario models, optimization functions,
  biological models, ...). Extend a hierarchy rather than adding free-standing scripts.
- Naming (enforced by MISS_HIT, see `.miss_hit`):
  - functions/scripts: `matRad_lowerCamelCase`; classes: `matRad_UpperCamelCase`;
    tests: `test_*`; test helpers: `helper_*`
  - variables, parameters, properties, methods: `lowerCamelCase` (all-caps acronyms allowed)
- Console output goes through the configuration object, never through bare
  `disp`/`fprintf`/`warning`/`error`:
  `matRad_cfg = MatRad_Config.instance(); matRad_cfg.dispInfo(...) / dispWarning(...) / dispError(...)`.
  Default parameters live in `MatRad_Config` too.
- Every `.m` file carries the matRad license header (copy it from any existing file and
  keep the year range current) and a docstring in the established format, since the
  API documentation is generated from it (see `docs/`).
- Keep functions small: `.miss_hit` limits cyclomatic complexity to 25, nesting depth to 6,
  parameters to 15 and file length to 1000 lines; line length is 150.
- Backwards compatibility matters for user-facing APIs (`pln`, `cst`, `ct`, `stf`, `dij`,
  `resultGUI` structs and public function signatures). If a signature has to change, keep
  accepting the old form with a deprecation warning where feasible.

## Pre-commit hooks and MISS_HIT

`.pre-commit-config.yaml` runs generic checks (whitespace, EOF, YAML, large files,
codespell, rstcheck for the docs) and the **MISS_HIT** tools on changed `.m` files:
`mh_style --fix` (formatting/naming), `mh_metric --ci` (complexity metrics) and `mh_lint`.
Install with `pip install pre-commit miss_hit` and `pre-commit install`.

Important: the hooks check the **whole file**, not just the diff. Many older files
predate MISS_HIT, so touching one line may produce a large amount of required (or
auto-applied) reformatting. Contributors may therefore **skip the hook**
(`git commit --no-verify` or `SKIP=mmh_style,mh_metric,mh_lint git commit ...`) to avoid
polluting a PR with unrelated churn — but the **portions you changed or added must
themselves satisfy the style, naming, metric and lint rules**. Run
`mh_style path/to/file.m`, `mh_metric --ci path/to/file.m` and `mh_lint path/to/file.m`
and fix everything that concerns your code. New files must pass cleanly in full.
Agents should never silently skip hooks; state it when doing so and why.

## Testing — test-driven development preferred

- Framework: [MOxUnit](https://github.com/MOxUnit/MOxUnit) (in `submodules/MOxUnit`),
  coverage via MOcov. Run everything with `matRad_runTests` from the repository root, or a
  subset with `matRad_runTests('test/optimization')`. Details in `test/README.md`.
- Write the test first (or alongside) for new functionality and for every bug fix; a bug
  fix without a regression test is incomplete. Place tests in the `test/` subfolder
  mirroring the source location, named `test_<feature>.m`, using the
  `test_functions = localfunctions(); initTestSuite;` header.
- Keep tests fast: use the small phantoms (`helper_createTestCt`, `matRad_PhantomBuilder`,
  `TG119`/`BOXPHANTOM` in `matRad/phantoms/`) and coarse resolutions; avoid GUI interaction
  and long dose calculations. Mark tests that require toolboxes/binaries so they skip
  cleanly (`moxunit_throw_test_skipped_exception`) instead of failing.
- Tests must pass on MATLAB R2022b **and** Octave. Be careful with exception identifiers
  in `assertExceptionThrown` on Octave.
- `matRad_cfg.setDefaultPropertiesForTesting()` is applied by the runner; do not rely on
  interactive defaults.
- Examples in `examples/` double as integration tests; if you change public behavior,
  check that the affected examples still run.

## Changelog, authorship and citation

- **Every user-relevant change is recorded in `CHANGELOG.md`** under `## [Unreleased]`
  in the appropriate `### Added` / `### Changed` / `### Deprecated` / `### Removed` /
  `### Fixed` / `### Security` section. The file follows
  [Keep a Changelog](https://keepachangelog.com/) and is maintained with the
  [kacl](https://github.com/mschmieder/python-kacl) conventions/tooling
  (`kacl-cli` can verify and add entries). Write entries from the user's perspective,
  name the affected functions/classes in backticks, and append them at the end of the
  respective list.
- **New contributors** add themselves to `AUTHORS.txt` (alphabetical by last name) and to
  the `authors` list in `CITATION.cff` (`family-names` / `given-names`, same order).
  Agents should not add names on their own; leave this to the human contributor.

## Documentation

- The documentation in `docs/` is built with **Sphinx**
  (`pip install -r docs/requirements.txt`, then `sphinx-build docs _doc` from the
  repository root; CI requires a warning-free build). Sections: `setup/`, `quickstart/`,
  `guide/` (technical guide), `datastructures/`, `overview/`, and `api/` (auto-generated
  from docstrings with `sphinxcontrib-matlabdomain`).
- Update the relevant `.rst` pages when behavior, options (`pln.prop*` fields), data
  structures or workflows change. New modules/classes must be reachable from `docs/api/`
  (classes in `@` folders need explicit `autoclass` entries, `automodule` does not descend
  into them).
- Docstrings are reStructuredText; malformed parameter tables render as prose, so build
  the docs (or at least run `rstcheck`) after touching them.
- Larger features should also get an example in `examples/` and a mention in the guide.

## Workflow checklist for a change

1. Understand the affected hierarchy/data structures; look for an existing base class or
   compat helper before writing new code.
2. Write or extend tests in `test/`; run `matRad_runTests('<folder>')` in MATLAB and, if
   available, in Octave.
3. Implement in `matRad/` following the conventions above (OOP, naming, `MatRad_Config`
   logging, absolute paths, R2022b/Octave compatible).
4. Run the MISS_HIT tools on the touched files; fix everything in your portions.
5. Update `CHANGELOG.md`, the docs in `docs/`, and (for larger features) `examples/`.
6. Do not commit generated artifacts (`coverage.*`, `testresults.xml`, `_doc/`, `build/`,
   `sphinxlog.txt`, local `.mex*` builds) or large data files.
7. Open the pull request against `dev` with a descriptive summary.

## Release workflow

Releases follow [Semantic Versioning](https://semver.org/) (see the versioning notes at the
end of `CHANGELOG.md`: major = breaking, minor = new features, patch = bug fixes). The
version is hard-coded in `matRad/matRad_version.m` (`matRadVer.name/major/minor/patch`) and
this is the single source of truth: `docs/conf.py` parses it for the Sphinx version, so do
not duplicate the number elsewhere; tags are named `vX.Y.Z` and the standalone packaging workflow (`.github/workflows/package.yml`)
builds automatically on pushes to `master`, `rc/**` and on tags. Only maintainers release;
agents may prepare the individual steps but must not push tags or delete branches without
explicit approval.

1. **Create the release candidate branch** from an up-to-date `dev`:
   `git checkout dev && git pull && git checkout -b rc/vX.Y.Z && git push -u origin rc/vX.Y.Z`
2. **Open a pull request** from `rc/vX.Y.Z` onto `master` (`gh pr create --base master`).
   Summarize the `[Unreleased]` changelog section in the PR description.
3. **Iterate on the release candidate** until the PR is approved: address review comments,
   fix CI failures (tests on R2022b / latest MATLAB / Octave, docs build, packaging) and close
   gaps in the patch coverage reported by the coverage workflow. Only fixes go onto
   `rc/vX.Y.Z`; new features wait for the next release. Cherry-pick or merge the fixes back
   to `dev` if they are needed there earlier.
4. **On approval, bump the version**: set `major`/`minor`/`patch` (and `name` for a new
   major/minor release) in `matRad/matRad_version.m`, move the `[Unreleased]` entries in
   `CHANGELOG.md` into a new `## [X.Y.Z] - YYYY-MM-DD` section including the compare links at
   the bottom (`kacl-cli release X.Y.Z` does this), update `version`/`date-released` in
   `CITATION.cff`, commit (`Bump version to vX.Y.Z`) and push. Wait for CI to go green again
   and for a maintainer to **merge the PR** with a **merge commit, never squash**: feature PRs
   into `dev` are already squashed, so the release PR must preserve those commits to keep the
   history of `master` and `dev` shared and the merge back into `dev` conflict-free.
5. **After the merge, tag and sync**:
   ```
   git checkout master && git pull
   git tag -a vX.Y.Z -m "matRad vX.Y.Z"
   git checkout dev && git pull && git merge master
   git push origin dev vX.Y.Z
   git push origin --delete rc/vX.Y.Z
   ```
   Then create the GitHub release from the tag (`gh release create vX.Y.Z --notes-from-tag`
   or with the changelog section as notes), attach the standalone installers produced by the
   packaging workflow, and check that the documentation build for the tag succeeded.

The tag name must match the version hard-coded in `matRad_version.m` exactly and sit on the
merge commit on `master`: `matRad_version` computes the revision number as the number of
commits since tag `v<major>.<minor>.<patch>` (revision `0` = release build). Hence the bump
happens before the tag, and a missing or mismatched tag leaves every build with an unknown
revision.

This sequence is a good candidate for a project-local command/skill (e.g.
`.claude/commands/matRad-release.md` or a `matRad-release` script) that walks through the
steps interactively and stops before every push/tag/delete for confirmation.


# renv 0.11.0

* Fixed an issue where `renv::install(..., type = "binary")` would
  still attempt to install packages from sources in some cases. (#461)
  
* `renv` now always writes `renv/.gitignore`, to ensure that the appropriate
  directories are ignored for projects which initialize `git` after `renv`
  itself is initialized. (#462)

* R Markdown documents with the `.Rmarkdown` file extension are now parsed for
  dependencies.

* Fixed an issue where setting the `external.libraries` configuration option
  would trigger a warning. (#452)

* Improved handling of unicode paths on Windows. (#451)

* `renv::snapshot(project = <path>)` now properly respects `.gitignore` /
  `.renvignore` files, even when that project has not yet been explicitly
  initialized yet. (#439)
  
* The default value of the `synchronized.check` option has been changed from
  TRUE to FALSE.

* Fixed an issue where packages downloaded from Bitbucket and GitLab did not
  record the associated commit hash.

* Fixed an issue where attempting to install packages from GitLab could fail
  to install the correct version of the package. (#436)

* `renv::snapshot()` now preserves records in a lockfile that are only
  available for a different operating system. This should make it easier
  to share lockfiles that make use of platform-specific packages. (#419)

* `renv` better handles files that are removed during an invocation to
  `renv::dependencies()`. (#429)

* The configuration option `install.staged` has been renamed to
  `install.transactional`, to better reflect its purpose. `install.staged`
  remains supported as a deprecated alias.

* Fixed an issue where `renv` could fail to parse non-ASCII content on Windows.
  (#421)

* `renv::update()` gains the `exclude` argument, useful in cases where one
  would like to update all packages used in a project, except for a small
  subset of excluded packages. (#425)

* `renv::update()` now respects the project `ignored.packages` setting. (#425)

* Fixed an issue where RSPM binary URL transformations could fail for
  Ubuntu Trusty. (#423)

* `renv` now records the `OS_type` reported in a package's `DESCRIPTION` file
  (if any), and ignores packages incompatible with the current operating
  system during restore. (#394)

# renv 0.10.0

* `renv::install()` gains the `type` argument, used to control whether `renv`
  should attempt to install packages from sources (`"source"`) or using
  binaries (`"binary"`).

* `renv` now knows how to find and activate Rtools40, for R 4.0.0 installations
  on Windows.

* The `RENV_PATHS_PREFIX` environment variable can now be used to prepend an
  optional path component to the project library and global cache paths.
  This is primarily useful for users who want to share the `renv` cache across
  multiple operating systems on Linux, but need to disambigutate these paths
  according to the operating system in use. See `?renv::paths` for more details.
  
* Fixed an issue where `renv::install()` could fail for packages from GitHub
  whose DESCRIPTION files contained Windows-style line endings. (#408)

* `renv::update()` now also checks and updates any Bioconductor packages
  used within a project. (#392)

* `renv` now properly parses negated entries within a `.gitignore`; e.g.
  `!script.R` will indicate that `renv` should include `script.R` when
  parsing dependencies. (#403)

* Fixed an issue where packages which had only binaries available on a
  package repository were not detected as being from a package repository.
  (#402)

* Fixed an issue where calls of the form `p_load(char = <vctr>)` caused a
  failure when enumerating dependencies. (#401)

* Fixed an issue where `renv::install()` could fail when multiple versions
  of a package are available from a single repository, but some versions of
  those packages are incompatible with the current version of R. (#252)

* Fixed an issue where downloads could fail when the associated pre-flight
  HEAD request failed as well. (#390)

* Fixed an issue where empty records within a DESCRIPTION file could cause
  `renv::dependencies()` to fail. (#382)

* renv will now download binaries of older packages from MRAN when possible.

* renv will now attempt to re-generate the system library sandbox if it is
  deleted while a session is active. (#361)

* Fixed an issue where Python packages referenced using `reticulate::import()`
  were incorrectly tagged as R package dependencies. Similarly, `renv` now only
  considers calls to `modules::import()` if those calls occur within a call to
  `modules::module()`. (#359)

* `renv::scaffold()` now also generates a lockfile when invoked. (#351)

* The argument `confirm` has been renamed to `prompt` in all places where it
  is used. `confirm` remains supported for backwards compatibility, but is no
  longer explicitly documented. (#347)

* The continuous integration `renv` vignette now also contains a template for
  using `renv` together with GitLab CI. (#348, @artemklevtsov)

* `renv` now properly resets the session library paths when calling
  `renv::deactivate()` from within RStudio. (#219)
  
* `renv::init()` now restores the associated project library when called in a
  project containing a lockfile but no project library nor any pre-existing
  project infrastructure.

* Fixed an issue on Windows where attempts to download packages from package
  repositories referenced with a `file://` scheme could fail.
  
* The configuration option `dependency.errors` has been added, controlling how
  errors are handled during dependency enumeration. This is used, for
  example, when enumerating dependencies during a call to `renv::snapshot()`.
  By default, errors are reported, and (for interactive sessions) the user is
  prompted to continue. (#342)
  
* `renv::dependencies()` gains two new arguments: the `progress` argument
  controls whether `renv` reports progress while enumerating dependencies,
  and `errors` controls how `renv` handles and reports errors encountered
  during dependency discovery. The `quiet` argument is now soft-deprecated,
  but continues to be supported for backwards compatibility. Specifying
  `quiet = TRUE` is equivalent to specifying `progress = FALSE` and
  `errors = "ignored"`. Please see the documentation in `?dependencies`
  for more details. (#342)
  
* The environment variable `RENV_PATHS_LIBRARY_ROOT` can now be set, to
  instruct `renv` to use a particular directory as a host for any project
  libraries that are used by `renv`. This can be useful for certain cases
  where it is cumbersome to include the project library within the project
  itself; for example, when developing an R package. (#345)

* The code used to bootstrap `renv` (that is, the code used to install `renv`
  into a project) has been overhauled. (#344)

* `renv` no longer unsets an error handler set within the user profile when
  loading a project. (#343)

* `renv` gains the "explicit" snapshot type, wherein only packages explicitly
  listed as dependencies within the project `DESCRIPTION` file (and those
  package's transitive dependencies) will enter the lockfile when
  `renv::snapshot()` is called. (#338)

* `renv` will now transform RSPM source URLs into binary URLs as appropriate,
  allowing `renv` to use RSPM's binary repositories during restore. See
  `?config` for more details. (#124)

* `renv` will now infer a dependency on `hexbin` in projects that make
  use of the `ggplot2::geom_hex()` function.

* `renv` now tries to place Rtools on the PATH when a package is installed
  with the `install.packages()` hook active. (#335)

# renv 0.9.3

* Fixed an issue where attempts to specify `RENV_PATHS_RTOOLS` would
  be ignored by `renv`. (#335)

* Fixed an issue where downloads could fail when using the `wininet`
  downloader, typically with a message of the form
  "InternetOpenUrl failed: 'The requested header was not found'".

* `renv` better handles projects containing special characters on Windows.
  (#334)

* `renv` better handles unnamed repositories. (#333)

* `renv` gains the config option `hydrate.libpaths`, allowing one to control
  the library paths used by default for `renv::hydrate()`. (#329)

* `renv::hydrate()` gains the `sources` argument, used to control the library
  paths used by `renv` when hydrating a project. (#329)
  
* `renv` now sandboxes the system library by default on Windows.

* `renv` now validates that the Xcode license has been accepted before
  attempting to install R packages from sources. (#296)
  
* The R option `renv.download.override` can now be used to override the
  machinery used by `renv` when downloading files. For example, setting
  `options(renv.download.override = utils::download.file)` would instruct
  `renv` to use R's own downloader when downloading files from the internet.
  This can be useful when configuration of `curl` is challenging or
  intractable in your environment, or you've already configured the base
  R downloader suitably.

* `renv::use_python("~/path/to/python")` now works as expected.

* `renv` now properly expands `R_LIBS_SITE` and `R_LIBS_USER` when set within a
  startup `.Renviron` file. (#318)

* The `renv.download.headers` option can now be used to provide arbitrary HTTP
  headers when downloading files. See the **Authentication** section in
  `vignette("renv")` for more details. (#307)

* `renv` gains the project setting `package.dependency.fields`, for controlling
  which fields in an R package's `DESCRIPTION` file are examined when
  discovering recursive package dependencies. This can be useful when you'd like
  to instruct `renv` to track, for example, the `Suggests` dependencies of the
  packages used in your project. (#315)

* `renv` now better handles repositories referenced using file URIs.

* Packages installed from GitHub using `renv::install()` will now also have
  `Github*` fields added, in addition to the default `Remote*` fields. This
  should help fix issues when attempting to deploy projects to RStudio Connect
  requiring packages installed by `renv`. (#397)
  
* `renv` now prefers using a RemoteType field (if any) when attempting to
  determine a package's source. (#306)

* `renv` gains a new function `renv::scaffold()`, for generating `renv` project
  infrastructure without explicitly loading the project. (#303)

* `renv` now updates its local `.gitignore` file, when part of a git repository
  whose git root lives in a parent directory. (#300)

# renv 0.9.2

* Fixed an issue in invoking `find` on Solaris.

# renv 0.9.1

* Fixed an issue in invoking `cp` on Solaris.

# renv 0.9.0

* `renv` gains a new function `renv::record()`, for recording new packages
  within an existing lockfile. This can be useful when one or more of the
  recorded packages need to be modified for some reason.

* An empty `.renvignore` no longer erroneously ignores all files within a
  directory. (#286)

* `renv` now warns if the version of `renv` loaded within a project does not
  match the version declared within the `renv` autoloader. (#285)

* `renv` gains a new function `renv::run()`, for running R scripts within
  a particular project's context inside an R subprocess. (#126)

* The algorithm used by `renv` for hashing packages has changed. Consider 
  using `renv::rehash()` to migrate packages from the old `renv` cache to
  the new `renv` cache.

* `renv::status()` now reports packages which are referenced in your project
  code, but are not currently installed. (#271)

* `renv` is now able to restore packages with a recorded URL remote. (#272)

* `renv::dependencies()` can now parse R package dependencies used as custom
  site generator in an Rmd yaml header. (#269, @cderv)

* `renv` now properly respects a downloader requested by the environment
  variable `RENV_DOWNLOAD_FILE_METHOD`.

* `renv` no longer sources the user profile (normally located at `~/.Rprofile`)
  by default. If you desire this behavior, you can opt-in by setting
  `RENV_CONFIG_USER_PROFILE = TRUE`; e.g. within your project or user
  `.Renviron` file. (#261)

* `renv::restore()` gains the `packages` argument, to be used to restore
  a subset of packages recorded within the lockfile. (#260)

* `renv` now tries harder to preserve the existing structure in infrastructure
  files (e.g. the project `.Rprofile`) that it modifies. (#259)

* `renv` now warns if any Bioconductor packages used in the project appear
  to be from a different Bioconductor release than the one currently active
  and stored in the lockfile. (#244)

* `renv` now normalizes any paths set in the `RENV_PATHS_*` family of
  environment variables when `renv` is loaded.

* Fixed an issue where `renv` would not properly clean up after a failed
  attempt to call `Sys.junction()`. (#251)

* Fixed an issue where `renv` would, in some cases, copy rather than link from
  the package cache when the library path had been customized with the
  `RENV_PATHS_LIBRARY` environment variable. (#245)

* The method `renv` uses when copying directories can now be customized. When
  copying directories, `renv` now by default uses `robocopy` on Windows, and
  `cp` on Unix. This should improve robustness when attempting to copy files
  in some contexts; e.g. when copying across network shares.

* `renv` now tracks the version of Bioconductor used within a project
  (if applicable), and uses that when retrieving the set of repositories
  to be used during `renv::restore()`.

* `renv::dependencies()` can now parse R package dependencies declared and
  used by the `modules` package. (#238, @labriola)

* Fixed an issue where `renv::restore()` could fail in Docker environments,
  usually with an error message like 'Invalid cross-device link'. (#243)

* `renv::install()` disables staged package install when running with the
  Windows Subsystem for Linux. (#239)

# renv 0.8.3

* `renv::dependencies()` gains a new argument `dev`, indicating whether
  development dependencies should also be included in the set of discovered
  package dependencies. By default, only runtime dependencies will be reported.

* `renv` has gained the function `renv::diagnostics()`, which can occasionally
  be useful in understanding and diagnosing `renv` (mis)behaviors.

* `renv::equip()` can now be used on macOS to install the R LLVM toolchain
  normally used when compiling packages from source. `renv` will also use
  this toolchain as appropriate when building packages from source.

* `renv::install()` now provides a custom Makevars when building packages on
  macOS with Apple Clang, to avoid issues due to the use of '-fopenmp' during
  compilation.

* `renv::install()` now respects explicit version requests when discovered
  in a project's DESCRIPTION file. (#233)

* Fixed an issue where `renv:::actions()` would fail to report any actions if
  the project lockfile was empty. (#232)

* When using `renv` for R package development, `renv` will no longer attempt to
  write the package being developed to the lockfile. (#231)

* Fixes for checks run on CRAN.

* renv will now search for Rtools in more locations. (#225)

* `renv::load()` now ensures that the version of `renv` associated with
  the loaded project is loaded when possible. In addition, experimental
  support for switching between projects with `renv::load()` has been
  implemented. (#229)

* `renv::dependencies()` no longer treats folders named with the extension
  `.Rmd` as though they were regular files. (#228)

* It is now possible to install source packages contained within `.zip`
  archives using `renv::install()`.

* Fixed an issue where attempts to call `renv::restore()` with the path to the
  lockfile explicitly provided would fail. (#227)

# renv 0.8.2

* Further fixes for checks run on CRAN.

# renv 0.8.1

* Fixes for checks run on CRAN.

# renv 0.8.0

* Initial CRAN release.

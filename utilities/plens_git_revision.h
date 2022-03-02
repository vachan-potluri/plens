#ifndef PLENS_GIT_REVISION_H
#define PLENS_GIT_REVISION_H

// This file will be used by `plens/CMakeLists.txt` (the main cmake script) to generate
// `plens_git_revision.h` in this directory. The latter file can be used to print git-related info
// from the code.

/**
 * Name of the local git branch of the source directory.
 */
#define PLENS_GIT_BRANCH "master"

/**
 * Full sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_REVISION "a8b7e2cb7e1cf0ae147c0a7104bda3ccd8864c58"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "a8b7e2c"

#endif

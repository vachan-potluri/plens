#ifndef PLENS_GIT_REVISION_H
#define PLENS_GIT_REVISION_H

// This file will be used by `plens/CMakeLists.txt` (the main cmake script) to generate
// `plens_git_revision.h` in this directory. The latter file can be used to print git-related info
// from the code.

/**
 * Name of the local git branch of the source directory.
 */
#define PLENS_GIT_BRANCH "subcell_interpolation"

/**
 * Full sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_REVISION "1a258973492b76b5b2f62265e0877d2eaa24ccf2"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define PLENS_GIT_SHORTREV "1a25897"

#endif

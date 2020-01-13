# HERA Validation Repository

**Archive of formal software pipeline validation tests**

This repository contains code and documentation for formal
tests of the HERA software pipeline. Tests are typically
performed and documented as Jupyter notebooks and are
archived to provide a long-standing account of the accuracy
of the pipeline as it evolves. Directory structures define the
broad kinds of tests performed.

## Mission Statement

The validation group seeks to validate the HERA data pipeline
software and algorithms by testing the specific software against
simulations where the expected output is well understood
theoretically.The group also helps to develop and define
increasingly sophisticated simulations on which to build an
end-to-end test and validation of the HERA pipeline.

## Structure of this repository

The validation effort seeks to verify the HERA software pipeline
through a number of well-defined steps of increasing complexity.
Each of these steps (called **major step**s or just **step**s in this
repository) reflects a broad validation concern or a specific 
element of the pipeline. For example, **step** 0 seeks to validate
just the ``hera_pspec`` software when given a well-known white-noise
P(k)-generated sky. 

Within each **step** exists the possibility of a set of variations 
(called **minor variation**s or just **variation**s in this repo). For 
example, variations for **step** 0 may be to generate flat-spectrum P(k)
and non-flat P(k). 

Finally, each combination of **step**-**variation** has the potential to incur
several staged tests or trials (we call them **trial**s in the repo). 

Importantly, failing **trial**s _will not be removed/overwritten_ in this
repo. Each formally-run trial is archived here for posterity. 

Thus the structure for this repo is as follows: Under the ``test-series``
directory, a number of directories labelled simply with their corresponding
**step** number are housed. Within each of these directories, each actual 
**trial** is presented as a notebook labelled ``test-<step>.<variation>.<trial>.ipynb``.

All **step**s, **variation**s and **trial**s are assigned increasing numerical
values. Generally, these values are increasing (from 0) in order of time/complexity.

In addition to the trial notebooks in these directories, each directory will
contain a ``README.md`` which lists the formal goals and conditions of each of
its **variation**s. 

Finally, each **variation** will be represented as a specific Github _project_,
in which the progress can be tracked and defined. Each project should receive 
a title which contains the **step**.**variation** identifier as well as a brief
description.

### Writing a validation test notebook

We have provided a template notebook which should serve as a starting
place for creating a validation notebook. The template is self-describing,
and has no intrinsic dependencies. All text in the notebook surrounded
by curly braces are meant to be replaced.

The template can be slightly improved/cleaned if you use jupyter notebook
extensions -- in particular the ToC and python-markdown extensions. The
first allows a cleaner way to build a table of contents (though one is
already included), and the latter allows using python variables in
markdown cells. This makes the writing of git hashes and versions simpler,
and means for example that the execution time/date can be written directly
into a markdown cell. 

## Project Plan
A semi-up-to-date version of the project plan is found at 
[project_table.md](https://github.com/HERA-Team/hera-validation/blob/project-table/project_table.md).

To create a simple tabulated version of the Project Plan, download the repo, save a
[personal access token](https://github.com/settings/tokens) to a file called `.pesonal-github-token`,
(ensure there is no trailing "\n" in the file)
and run `make_project_table.py` at the root directory. 
Note that you will need python 3.4+ and the `pygithub` code to run this script (`pip install pygithub`).

Do this on the `project-table` branch to update our "current" table linked above.

## Project Organization
To track where the validation effort is currently at, we rely on a particular structure
in this repository. We first ask you to review the overall structure of the test series
(above)[#structure-of-this-repository] as well as in the 
[project table](https://github.com/HERA-Team/hera-validation/blob/project-table/project_table.md).
This will indicate what the relevant step number is for any proposed test. 

### Proposing a new test
To propose a new test, create a 
(new test-proposal issue)[https://github.com/HERA-Team/hera-validation/issues/new?assignees=&labels=formal-test&template=test_proposal.md&title=Step+X.X%3A+%3CShort+Title%3E].
Fill out the information in the issue template with as much detail as possible, in 
particular, which simulation components will be required (and their parameters), and
which pipeline components will be required. If you have access, you can also update the
tags applied to the issue. In particular, there are specific tags for the following:
* Simulator used
* Simulated components (eg. GLEAM, GSM, gains, RFI, ...)
* Pipeline components required (eg. pspec, abscal, ...)
It is these tags (not the description in the actual issue body) that are displayed in our
project table summary, so this is important. However, we will add these tags during our
weekly review if you are unable to do it.

Finally, you should (if you can) apply the `status:proposed` tag, to inform the team that
this is a new test to be considered for formal inclusion in the test series. Once reviewed
at our weekly meeting, we will either promote this to `status:accepted` or demote it
to `status:rejected`. 

If you are willing to do the work to perform the test you have proposed, you may also
assign yourself. You may also contact any others you think might be able to help you
and assign them. 

### Weekly Review
As has already been mentioned, during the weekly validation meeting, we will assess any
new issues with the tag `status:proposed`, and update the status tag, assignees and 
metadata tags (simulators, components, pipelines, etc.). It will be particularly 
helpful if the proposer of the test could be present at the weekly meeting directly 
following the proposal.  

In addition to updating these tags, the team will decide on two things:
1. Which (project)[projects] to formally add the proposed test to (these correspond to 
   the various validation steps)
2. Which (milestone)[milestones] to add the proposed test to (if any). These represent 
   broader goals of the project. For example, the first milestone is 
   (H1C_IDR2)[https://github.com/HERA-Team/hera-validation/milestone/4]. This sets
   basic priorities for the team.
   
### Starting/Running a Test
To start working on a test, first create a branch (preferrably named after the test
step explicitly), and proceed to edit a (validation template)[notebook-template.ipynb].
After your first commit/push, create a _draft_ PR, which in the body explicitly 
references the original issue proposing the test. From this time, the original
test proposal will be considered frozen, and all development of the ideas, and discussion,
will be had in the PR itself. It is important that the original issue be linked explicitly
with a line such as "Fixes #37", in order for the issue tracking to work as we want. 

In detail, doing so will automatically set the status tag on the original issue to
`status:active`, and also automatically move it from "To Do" to "In Progress" in the 
relevant project. 

Do *not* set the Project/Milestone for the PR (but do assign anyone who is part of the
execution of the PR). Doing so duplicates the step in the Project tracker.

While your PR is in draft form, we will attempt to discuss updates on it weekly at our
meeting (if time permits). When it is ready for final review, set it to a non-draft PR,
and we will review it _as a team_ at the next available meeting. We will also assign
a single final reviewer who will give the thumbs up for merging.

### Tracking Progress
If you wish to track the progress of the validation effort, the most in-depth information
can be determined from the (projects)[projects] page. Each Project is a major validation
"Step", and opening any particular project will show which sub-steps are to do, in progress
or completed. 

Do remember that projects are irrespective of overall goal/milestone, and so you should
filter them by the current milestone as well.  
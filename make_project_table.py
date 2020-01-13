"""
Make the Project Plan table by accessing GitHub issue/PR information.

Requires python 3.4+ and pygithub to be installed.
"""
from github import Github
import re
import os
import sys

table_header = """Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
"""

if __name__=="__main__":
    if len(sys.argv) > 1:
        GH_TOKEN = sys.argv[-1]
    else:
        if not os.path.exists('.personal-github-token'):
            raise ValueError(
                "You do not yet have a personal access token for github installed. Create one at "
                "https://github.com/settings/tokens "
                "and paste it into the file .personal-github-token (notice leading .) in this "
                "directory.")

        with open(".personal-github-token") as fl:
            GH_TOKEN = fl.read()

    g = Github(GH_TOKEN)
    repo = g.get_repo("HERA-Team/hera-validation")

    projects = repo.get_projects()
    step_number_pattern = re.compile(r"(-?\d+(\.\d+){0,2})")

    tables = """## Project Plan
Here follows a formal plan for envisaged tests, automatically generated from our Projects. 
This list is *not final*: it may be extended or 
modified at any time. However, it provides a reasonable snapshot of the current status and plans of 
the validation team.

What the symbols mean:

Sym. | Meaning | Sym. | Meaning | Sym.  |Meaning
-------| ----- | ----- | ---- | ----- | ------
:egg:      | Idea in gestation  | :hammer:   | Currently being worked on | :thinking: | PR submitted... 
:heavy_check_mark: | Passed     | :x:        | Failed | :watch:    | Not required for H1C IDR2

"""

    for proj in sorted(projects, key=lambda p: p.name):

        cols = proj.get_columns()

        steps = sum([list(col.get_cards()) for col in cols], [])
        steps = [step.get_content() for step in steps]

        step_dict = {}

        for step in steps:
            if step.__class__.__name__ not in ["Issue", "PullRequest"]:
                print("card ", step, " has class", step.__class__.__name__)
                continue
            # get status
            if step.milestone is None or step.milestone.title != "H1C IDR2":
                status = ":egg:"
                status_num = 0
            elif step.__class__.__name__ == "Issue" and step.state == 'open':
                status = ":hammer:"
                status_num = 1
            elif step.__class__.__name__ == "PullRequest" and step.state == "open":
                status = ":thinking:"
                status_num = 2
            elif step.state == "closed":
                status = ":heavy_check_mark:"
                status_num = 3
            else:
                raise ValueError("status not recognized for ", step)

            # Get step number
            step_number = step_number_pattern.search(step.title).group()
            major_step = ".".join(step_number.split(".")[:2])  # only up to X.X
            description = step.title.split(step_number)[-1].split(":")[-1].split('--')[-1]
            if description:
                description = ": " + description

            # Get title
            title = f"[{major_step}]({step.url})"

            # Get simulators used
            labels = step.labels
            sims = ", ".join([lbl.name.split(":")[-1] for lbl in labels if lbl.name.split(":")[0] == 'simulator'])

            simcmp = ", ".join([lbl.name.split(":")[-1] for lbl in labels if lbl.name.split(":")[0] == 'simcmp'])
            anlcmp = ", ".join([lbl.name.split(":")[-1] for lbl in labels if lbl.name.split(":")[0] == 'pipeline'])
            assigned = ", ".join([f"[@{assgn.login}]({assgn.url})" for assgn in step.assignees])

            if major_step in step_dict:
                if status_num <= step_dict[major_step]['status_num']:
                    continue

            step_dict[major_step] = dict(
                status = status,
                title = title,
                description=description,
                sims = sims,
                simcmp = simcmp,
                anlcmp = anlcmp,
                assigned = assigned,
                status_num = status_num
            )

        tables += f"### [{proj.name}]({proj.url})\n"
        tables += f"{proj.body}\n\n"
        tables += table_header
        for major_step in sorted(step_dict.keys()):
            d = step_dict[major_step]
            tables += f"{d['status']}  |  {d['title']}{d['description']}  | {d['sims']}  | {d['simcmp']} | {d['anlcmp']} | {d['assigned']}  |\n"

        tables += "\n\n"

    with open("project_table.md", 'w') as fl:
        fl.write(tables)



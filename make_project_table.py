import sys
from github import Github
import re

GH_TOKEN = "270073b246bb54cc9c322c07368c57131fdec866"

table_header = """Status     | #    | Simulator(s) | Sim. Components | Analysis Components | Assigned |
-----------| -----|--------------|-----------------|---------------------|----------|
"""

if __name__=="__main__":
    g = Github(GH_TOKEN)
    repo = g.get_repo("HERA-Team/hera-validation")

    projects = repo.get_projects()

    # issues = repo.get_issues(state='all')
    # prs = repo.get_pulls(state='all')
    #
    # # Isolate issues/prs that define formal tests/steps
    # steps = [issue for issue in everything if "formal-test" in [lbl.name for lbl in issue.labels]]
    #
    # everything = [issue for issue in issues] + [pr for pr in prs]
    # This finds anything with the pattern X.X.X in a string (i.e. the STEP.MAJOR part of a
    # step identifier)
    step_number_pattern = re.compile(r"(-?\d+(\.\d+){0,2})")

    tables = ""
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

        tables += f"[{proj.name}]({proj.url}): {proj.body}\n"
        tables += table_header
        for major_step in sorted(step_dict.keys()):
            d = step_dict[major_step]
            tables += f"{d['status']}  |  {d['title']}{d['description']}  | {d['sims']}  | {d['simcmp']} | {d['anlcmp']} | {d['assigned']}  |\n"

        tables += "\n\n"

    with open("project_table.md", 'w') as fl:
        fl.write(tables)



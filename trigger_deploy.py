import subprocess
import py


about = {}
with open("o3seespy/__about__.py") as fp:
    exec(fp.read(), about)


version = about['__version__']

failures = py.test.cmdline.main(['tests/'])
if failures == 0:
    subprocess.check_call(["git", "tag", version, "-m", "version %s" % version])
    subprocess.check_call(["git", "push", "--tags"])

# git push --tags origin pypi
# git tag 0.5.2 -m "version 0.5.2"
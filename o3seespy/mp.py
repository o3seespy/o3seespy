from o3seespy.base_model import opy


def barrier():
    opy.barrier()


def get_np():
    return opy.getNP()


def get_pid():
    return opy.getPID()


def partition():
    return opy.partition()

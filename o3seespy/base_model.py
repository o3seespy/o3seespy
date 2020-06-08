from collections import OrderedDict
try:
    import custom_openseespy.opensees as opy
except ModuleNotFoundError:
    import openseespy.opensees as opy
from o3seespy import exceptions
from o3seespy import extensions


class OpenSeesObject(object):
    op_base_type = "<not-set>"  # used to call opensees module
    op_type = "<not-set>"  # string name given to object in opensees module
    _tag = None
    _name = None
    _parameters = None

    @property  # deliberately no setter method
    def tag(self):
        return self._tag

    @property
    def name(self):
        return self._name

    def to_process(self, osi):
        if osi.state == 0:
            self.to_opensees()
        if osi.state == 1:
            osi.to_commands(self.to_commands())
        elif osi.state == 2:
            osi.to_dict(self)
        elif osi.state == 3:
            osi.to_commands(self.to_commands())
            self.to_opensees()
        elif osi.state == 4:
            osi.to_commands(self.to_commands())

    def to_opensees(self):
        try:
            try:
                return getattr(opy, self.op_base_type)(*self.parameters)
            except opy.OpenSeesError as e:
                com = extensions.to_commands(self.op_base_type, self.parameters)
                raise ValueError('{0} caused error "{1}"'.format(com, e))

        except SystemError as e:
            if None in self.parameters:
                print(self.parameters)
                raise exceptions.ModelError("%s of type: %s contains 'None'" % (self.op_base_type, self.op_type))
            else:
                raise SystemError(e)
        except AttributeError as e:
            print(e)
            print('opensees.{0}({1}) caused error "{2}"'.format(self.op_base_type,
                                                                               ','.join(str(x) for x in self.parameters),
                                                                                        e))
            raise exceptions.ModelError("op_base_type: '%s' does not exist in opensees module" % self.op_base_type)

    def to_commands(self):
        return extensions.to_commands(self.op_base_type, self.parameters)

    @property
    def parameters(self):
        return self._parameters

    @property
    def type(self):
        return self.op_type

    @property
    def base_type(self):
        return self.op_base_type

    def to_dict(self, export_none=False):
        outputs = OrderedDict()
        for item in self.__dict__:
            if '_' == item[0]:  # do not export private variables
                continue
            value = self.__getattribute__(item)
            if not export_none and value is None:
                continue
            outputs[item] = collect_serial_value(value)
        return outputs


class OpenSeesMultiCallObject(OpenSeesObject):
    _multi_parameters = None

    @property
    def parameters(self):
        return self._multi_parameters[-1]

    @property
    def multi_parameters(self):
        return self._multi_parameters


def collect_serial_value(value):
    if isinstance(value, str):
        return value
    elif isinstance(value, int):
        return value
    elif hasattr(value, "to_dict"):
        return value.to_dict()
    elif hasattr(value, "__len__"):
        tolist = getattr(value, "tolist", None)
        if callable(tolist):
            value = value.tolist()
            return value
        else:
            if hasattr(value, "tag"):
                return value.tag
            elif hasattr(value, "to_dict"):
                value = value.to_dict()
                return value
            else:
                values = []
                for item in value:
                    values.append(collect_serial_value(item))
                return values
    else:
        return value



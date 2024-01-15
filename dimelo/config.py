import platform
from pathlib import Path
from dataclasses import dataclass, field, InitVar


@dataclass
class DependencyConfig:
    """Config definitions for depenency locations; only need to pass the *_dir args to the constructor."""
    modkit_dir: InitVar[Path | str]
    htslib_dir: InitVar[Path | str]

    modkit_exe: Path = field(init=False)
    bgzip_exe: Path = field(init=False)
    tabix_exe: Path = field(init=False)

    def __post_init__(self,
                      modkit_dir,
                      htslib_dir):
        modkit_dir = Path(modkit_dir)
        htslib_dir = Path(htslib_dir)
        self.modkit_exe = modkit_dir / 'modkit'
        self.bgzip_exe = htslib_dir / 'bgzip'
        self.tabix_exe = htslib_dir / 'tabix'

# Mapping of {platform.node(): DependencyConfig}
# THIS IS WHAT TO ACTUALLY UPDATE/MODIFY FOR NEW SYSTEMS
__dependency_configs = {
    'The-Kugel-MacBook.local': DependencyConfig(modkit_dir='/Users/jeremy/.cargo/bin',
                                                htslib_dir='/opt/homebrew/bin'),
    # 'Oberons-MacBook-Pro.local': DependencyConfig(modkit_dir='/Users/oberondixon-luinenburg/Documents/GitHub/modkit/target/release',
    #                                               htslib_dir='/opt/homebrew/bin'),
}


# Code-facing configuration setting
if platform.node() in __dependency_configs:
    print(f'Node {platform.node()} has a pre-defined set of executable directories')
    EXE_CONFIG = __dependency_configs[platform.node()]
elif platform.system()=='Linux':
    print('Default Linux configuration: ../dependencies/linux for executables.')
    EXE_CONFIG = DependencyConfig(modkit_dir='../dependencies/linux/modkit',
                                  htslib_dir='../dependencies/linux/htslib/bin')
elif platform.system()=='Darwin':
    print('Default Mac OSX configuration: ../dependencies/macosx for executables.')
    EXE_CONFIG = DependencyConfig(modkit_dir='../dependencies/macosx/modkit',
                                  htslib_dir='../dependencies/macosx/htslib')
else:
    raise(f'System {platform.system()} and node {platform.node()} not found in config.py' )

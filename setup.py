from setuptools import setup, find_packages


with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='flupre',
    version='1.0.0',
    author='viralInformatics',
    author_email='pys2013@hnu.edu.cn',
    description='A tool for predicting influenza virus hosts and virulence level from sequence data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/viralInformatics/FluRisk',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'flupre=flupre.flupre:main',
            'barplot=flupre.plot_bar:main',
            # 'radarplot=flupre.risk_assessment:main',
            'risk2tab=flupre.risk_assessment:main',
        ],
    },
    # package_data={
    #     '': ['data/*.*', 'model/*.*','tests/*.*', '18Mid/*.*', 'temp/*.*','app/*.*', 'querySeq/*.*','result/*.*'],
    # },
    # 所有想保存的非脚本的包必须是有init文件的package才行，没有init就不是package，那么这个package_data参数自然无用
    package_data={
        '': ['data/*', 'model/*','test/*', '18Mid/*', 'temp/*','app/*', 'querySeq/*','result/*'],
    },
    # package_data={
    #     '': ['flupre/*', 'fluhp/*','test/*', 'fluvp/*', 'convert_site/*'],
    # },
    install_requires=required,
    zip_safe=False,
 
    classifiers=[
        'Programming Language :: Python :: 3.6',
        # 'Programming Language :: Python :: 3.7',
        # 'Programming Language :: Python :: 3.8',
        # 'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent',
    ],
    # python_requires='>=3.6',
    python_requires='>=3.6, <3.8',
    keywords='influenza risk assessment tool',
)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
from datetime import datetime


def main_logo():
    print("""
                   ___   _      ___  __    ___ 
     _ __  _   _  / _ \ /_\    / __\/ _\  / _ \\
    | '_ \| | | |/ /_\///_\\\\  / /   \ \  / /_)/
    | |_) | |_| / /_\\\\/  _  \/ /___ _\ \/ ___/ 
    | .__/ \__, \____/\_/ \_/\____/ \__/\/     
    |_|    |___/                               

    Author: Jianghai Wang@BUAA
    Email: wang_jianghai@buaa.edu.cn

    Welcome to the pyGACSP program.
    """)


def pop_initialization_logo():
    print('''
      _____       _ _   _       _ _          _   _             
      \_   \_ __ (_) |_(_) __ _| (_)______ _| |_(_) ___  _ __  
       / /\/ '_ \| | __| |/ _` | | |_  / _` | __| |/ _ \| '_ \ 
    /\/ /_ | | | | | |_| | (_| | | |/ / (_| | |_| | (_) | | | |
    \____/ |_| |_|_|\__|_|\__,_|_|_/___\__,_|\__|_|\___/|_| |_|
                                                           

        Author: Jianghai Wang@BUAA   2024

        Initializing the population...
        
        Time: {}
        '''.format(datetime.now()))


def pop_iteration_logo():
    print('''
      _____ _                 _   _                  __ _             _   
      \_   \ |_ ___ _ __ __ _| |_(_) ___  _ __      / _\ |_ __ _ _ __| |_ 
       / /\/ __/ _ \ '__/ _` | __| |/ _ \| '_ \     \ \| __/ _` | '__| __|
    /\/ /_ | ||  __/ | | (_| | |_| | (_) | | | |    _\ \ || (_| | |  | |_ 
    \____/  \__\___|_|  \__,_|\__|_|\___/|_| |_|    \__/\__\__,_|_|   \__|
    
    
        Start population iteration...
        
        Time: {}
    '''.format(datetime.now()))


def finish_logo(gen_num):
    print('''
       __    __           _         ___                        _  _  _ 
      / / /\ \ \___  _ __| | __    /   \___  _ __   ___       / \/ \/ \\
      \ \/  \/ / _ \| '__| |/ /   / /\ / _ \| '_ \ / _ \     /  /  /  /
       \  /\  / (_) | |  |   <   / /_// (_) | | | |  __/    /\_/\_/\_/ 
        \/  \/ \___/|_|  |_|\_\ /___,' \___/|_| |_|\___|    \/ \/ \/   
                                                                 
        Work Done!
        
        Total generations: {}
        
        Time: {}
    '''.format(gen_num, datetime.now()))


def restart_logo(gen):
    print('''
       __           _             _          _  _  _ 
      /__\ ___  ___| |_ __ _ _ __| |_       / \/ \/ \\
     / \/// _ \/ __| __/ _` | '__| __|     /  /  /  /
    / _  \  __/\__ \ || (_| | |  | |_     /\_/\_/\_/ 
    \/ \_/\___||___/\__\__,_|_|   \__|    \/ \/ \/   
    
        RESTART from generation: {}
        
        Time: {}
    '''.format(gen, datetime.now()))


# def create_log():
#     import logging
#
#     logger = logging.getLogger()
#     logger.setLevel(logging.DEBUG)
#
#     file_handler = logging.FileHandler('running.log')
#
#     formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#     file_handler.setFormatter(formatter)
#
#     logger.addHandler(file_handler)


if __name__ == '__main__':
    pop_initialization_logo()
